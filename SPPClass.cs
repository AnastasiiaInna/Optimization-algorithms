using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Threading;
using System.Reflection;
using MathsLib;

namespace SPP
{
    class MyRandom : Random
    {
        private Random r;
        public MyRandom() : base((int)DateTime.Now.Ticks) { }
        /*public override int Next(int Seed)
        {
            return base.Next(Seed);
        }*/
    }

    class SPPClass
    {
        #region Types
        delegate void VoidFunc();
        delegate void VoidFunc<T1>(T1 p1);
        delegate void VoidFunc<T1, T2>(T1 p1, T2 p2);
        delegate void VoidFunc<T1, T2, T3>(T1 p1, T2 p2, T3 p3);
        public enum AlgType
        {
            atGeneticAlgorithm,
            atHillClimbing,
            atTabuSearch,
            atScatterSearch
        }

        public enum MarkType
        {
            mtNew,
            mtOld
        }
        #region MemoryPopulationElement
        public class MemoryPopulationElement : System.Collections.IEnumerable
        {
            #region MemoryPopulationElement Fields
            private double _Fitness;
            private double _fitness;
            private double _unfitness;
            public byte[] Value;
            private int _Count;
            private MemoryPopulationClass _parent;
            private int[] _RowsCoverage;
            private List<int> _CoveredRows;
            // private MarkType _mark = MarkType.mtNew;
            //private int[] _RowsCoverage;
            #endregion

            #region MemoryPopulationElement Properties
            public int[] RowsCoverage
            {
                get { return _RowsCoverage; }
                set { _RowsCoverage = value; }
            }

            public List<int> CoveredRows
            {
                get { return _CoveredRows; }
                set { _CoveredRows = value; }
            }

            public MemoryPopulationClass Parent
            {
                get { return _parent; }
                set { _parent = value; }
            }

            public double Fitness
            {
                get { return _Fitness; }
                set { _Fitness = value; }
            }

            public double fitness
            {
                get { return _fitness; }
                set { _fitness = value; }
            }

            public double unfitness
            {
                get { return _unfitness; }
                set { _unfitness = value; }
            }

            public byte[] x
            {
                get { return Value; }
                set { Value = value; }
            }

            public int Count
            {
                get { return _Count; }
                set
                {
                    if (value >= 0)
                    {
                        _Count = value;
                        FillRandom();
                        Process();
                    }
                }
            }

            public byte this[int Index]
            {
                get { return Value[Index]; }
                set { Value[Index] = value; }
            }

            public int Index
            {
                get;
                set;
            }

            // public MarkType Mark
            //{
            //    get { return _mark; }
            //    set { _mark = value; }
            //}
            #endregion

            #region MemoryPopulationElement Constructors
            public MemoryPopulationElement(MemoryPopulationClass Parent, int Count, bool ToProcess = true)
            {
                this._parent = Parent;
                _Fitness = 0;
                _Count = Count;
                FillRandom();
                if (ToProcess) Process();
            }
            public MemoryPopulationElement(MemoryPopulationClass Parent, byte[] arr, int Count = 0, bool ToProcess = true)
            {
                this._parent = Parent;
                _Fitness = 0;
                _Count = Count;
                Fill(arr);
                if (ToProcess) Process();
            }

            public MemoryPopulationElement(MemoryPopulationClass Parent, byte[] x, bool ToProcess = true)
            {
                this._parent = Parent;
                _Fitness = 0; _Count = x.Length;
                Value = new byte[x.Length];
                x.CopyTo(Value, 0);
                if (ToProcess) Process();
            }
            public MemoryPopulationElement(MemoryPopulationElement p, bool FullCopy = false)
            {
                this._parent = p.Parent;
                this._Fitness = p.Fitness;
                this._unfitness = p.unfitness;
                this._fitness = p.fitness;
                this._Count = p.Count;
                //this._mark = p.Mark;
                this.Value = new byte[p.Count];

                p.x.CopyTo(this.Value, 0);

                if (FullCopy)
                {
                    this._RowsCoverage = new int[Parent.Parent.RowCount];
                    p.RowsCoverage.CopyTo(this._RowsCoverage, 0);

                    this._CoveredRows = new List<int>();
                    this._CoveredRows.AddRange(p.CoveredRows);
                }
            }
            #endregion

            #region MemoryPopulationElement Methods
            public MemoryPopulationElement()
            {
            }

            public MemoryPopulationElement Copy(bool FullCopy = false)
            {
                return new MemoryPopulationElement(this, FullCopy);
            }

            public void Process(bool GetFitnessOnly = false)
            {
                GetFitness();
                if (!GetFitnessOnly)
                {
                    GetRowsCoverage();
                    _CoveredRows = RowsCoverd();
                }
            }

            public void FillRandom()
            {
                MyRandom R = new MyRandom();
                Value = new byte[Count];
                for (int i = 0; i < _Count; i++)
                {
                    Value[i] = (byte)R.Next(2);
                }
            }

            public void Fill(byte[] arr)
            {
                //   MyRandom R = new MyRandom();
                Value = new byte[Count];
                for (int i = 0; i < _Count; i++)
                {
                    Value[i] = arr[i];
                }
            }

            private double GetFitness()
            {
                return _Fitness = SPPClass.Fitness(Value, Parent.Parent.A, Parent.Parent.c, out _fitness, out _unfitness);
            }

            public System.Collections.IEnumerator GetEnumerator()
            {
                return Value.GetEnumerator();
            }

            private void GetRowsCoverage()
            {
                _RowsCoverage = new int[Parent.Parent.RowCount];

                for (int i = 0; i < Parent.Parent.RowCount; i++)
                {
                    _RowsCoverage[i] = RowCovering(i);
                }
            }

            public int RowCovering(int RowIndex)
            {
                int res = 0;
                for (int i = 0; i < Count; i++)
                {
                    res += Parent.Parent.A[RowIndex][i] * this[i];
                }

                return res;
            }

            public List<int> RowsCoverd()
            {
                List<int> L = (from i in Parent.Parent.RowIndexes where RowsCoverage[i] > 0 select i).ToList();
                return L;
            }

            public double GetDistance(MemoryPopulationElement Other)
            {
                if (this._Count != Other.Count)
                {
                    throw new Exception("Individuals should have identical sizes!");
                }
                double Result = 0;

                for (int i = 0; i < _Count; i++)
                {
                    Result += (Value[i] - Other[i]) * (Value[i] - Other[i]);
                }
                Result = Math.Sqrt(Result);

                return Result;
            }

            public double GetDistance(MemoryPopulationClass MemoryPopulation)
            {
                return (from p in MemoryPopulation select this.GetDistance(p)).OrderBy(d => d).ToList()[0];
            }

            internal void SetParent(MemoryPopulationClass Parent, bool ToProcess)
            {
                _parent = Parent;
                if (ToProcess)
                {
                    Process(true);
                }
            }
            #endregion
        }
        #endregion
        #region MemoryPopulationClass
        public class MemoryPopulationClass : IEnumerable<MemoryPopulationElement>
        {
            #region MemoryPopulationClass Fields
            private List<MemoryPopulationElement> _memorypopulation = new List<MemoryPopulationElement>();

            private int _Count;
            private int _ColCount;

            private SPPClass _parent;
            #endregion

            #region PopulationClass Properties
            public SPPClass Parent
            {
                get { return _parent; }
                set { _parent = value; }
            }

            public List<MemoryPopulationElement> MemoryPopulation
            {
                get { return _memorypopulation; }
            }

            public int Count
            {
                get { return _Count; }
            }

            public int ColCount
            {
                get { return _ColCount; }
            }

            public MemoryPopulationElement this[int Index]
            {
                get { return _memorypopulation[Index]; }
                set { _memorypopulation[Index] = value; }
            }

            protected MemoryPopulationElement Last
            {
                get
                {
                    return this[_Count - 1];
                }

                set
                {
                    this[_Count - 1] = value;
                }
            }

            #endregion

            public MemoryPopulationClass(SPPClass Parent, int Count = 0)
            {
                this._parent = Parent;
                _ColCount = Parent.ColCount;
                _Count = Count;
                Generate();
            }

            public void Generate(byte[] p)
            {
                MyRandom R = new MyRandom();
                //for (int k = 0; k < Ms; k++)
                //{
                //    int ind;
                //    //  do
                //    // {
                //    ind = R.Next(E.Count);
                //    //ind = 
                //    //   } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

                //    ChangeIndex.Add(ind);
                //    E[ind] = (byte)(1 - E[ind]);//(byte)R.Next(2);
                //}

                _memorypopulation.Clear();// = new PopulationElement[Count];
                byte[] temp = new byte[p.Count()];
                temp = p.ToArray();
                //int k = R.Next(0, temp.Count() - 1);
                int k = (byte)R.Next(temp.Count() - 1);
                for (int i = 0; i < Count; i++)
                {
                    //k++;
                    //if (k == temp.Count() -1)
                    //    k--;
                    //else k ++;
                    //k = (byte)R.Next(temp.Count());    
                    //p[k] = (byte)Math.Abs(temp[k] - 1);
                    MemoryPopulationElement tmp;
                    do
                    {
                        {
                            p[k] = (byte)Math.Abs(temp[k] - 1);
                            tmp = new MemoryPopulationElement(this, p, ColCount);
                            if (Parent.UsingHeuristicOperator)
                                MemoryPopulationClass.HeuristicFeasibilityOperator(ref tmp);
                            k = (byte)R.Next(temp.Count() - 1);
                            //if (k == p.Count() - 1)
                            //    k--;
                            // else k++;
                            // if (i < Count / 2)
                            //     k++;
                            // else k -= (i - 1);
                        }
                    }
                    while (this.Exists(tmp));

                    _memorypopulation.Add(tmp);
                    _memorypopulation[i].Index = i;
                }
            }
            public MemoryPopulationClass(SPPClass Parent, int eps, byte[] p)
            {
                this._parent = Parent;
                _ColCount = Parent.ColCount;
                _Count = eps;
                Generate(p);
            }

            #region MemoryPopulationClass Methods
            public MemoryPopulationClass Add(MemoryPopulationElement p)
            {
                if (this.Exists(p)) return this;

                _memorypopulation.Add(p);
                MemoryPopulationElement tmp = _memorypopulation[_Count];// = p.Copy();
                //if (tmp.Parent.Parent == this.Parent)
                //    tmp.SetParent(this, false);
                //else 
                    tmp.SetParent(this, true);
                _Count++;

                return this;
            }

            public MemoryPopulationClass Add(MemoryPopulationClass OtherPopulation)
            {
                foreach (MemoryPopulationElement p in OtherPopulation)
                {
                    this.Add(p);
                }

                return this;
            }

            public MemoryPopulationClass Delete(MemoryPopulationElement p)
            {
                _memorypopulation.Remove(p);// = (from pp in _population where pp != p select pp).ToList();
                _Count--;
                return this;
            }

            public void Generate()
            {
                _memorypopulation.Clear();// = new PopulationElement[Count];

                for (int i = 0; i < Count; i++)
                {
                    MemoryPopulationElement tmp;
                    do
                    {
                        tmp = new MemoryPopulationElement(this, ColCount);
                        if (Parent.UsingHeuristicOperator)
                            MemoryPopulationClass.HeuristicFeasibilityOperator(ref tmp);
                    }
                    while (this.Exists(tmp));

                    _memorypopulation.Add(tmp);
                    _memorypopulation[i].Index = i;
                }
            }

            public double Fitness(int Index)
            {
                if (Index >= 0 && Index < Count)
                    return MemoryPopulation[Index].Fitness;
                else return Double.MaxValue;
            }

            public void Sort()
            {
                _memorypopulation = _memorypopulation.OrderBy(element => element.Fitness).ToList();
            }

            public static MemoryPopulationElement Crossover(MemoryPopulationElement P1, MemoryPopulationElement P2)
            {
                if (P1.Count == P2.Count && P1.Parent == P2.Parent)
                {
                    int count = P1.Count;
                    MemoryPopulationElement Mask = new MemoryPopulationElement(P1.Parent, count);
                    MemoryPopulationElement Child = new MemoryPopulationElement(P1.Parent, count, false);

                    for (int i = 0; i < count; i++)
                    {
                        Child[i] = (Mask[i] == 0) ? (P1[i]) : (P2[i]);
                    }

                    Child.Process();
                    return Child;
                }
                else return null;
            }


            public static void Mutation(ref MemoryPopulationElement E, bool DynamicMutation = false)
            {
                //     int[] Theta = new int[E.Parent.Parent.RowCount];

                int i, p;
                double eps = 0.5;
                int Md = 5, Ms = 3;
                List<int> ChangeIndex = new List<int>();

                //============Static Mutation================
                MyRandom R = new MyRandom();
                for (int k = 0; k < Ms; k++)
                {
                    int ind;
                    //  do
                    // {
                    ind = R.Next(E.Count);
                    //ind = 
                    //   } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

                    ChangeIndex.Add(ind);
                    E[ind] = (byte)(1 - E[ind]);//(byte)R.Next(2);
                }

                E.Process();

                //if (DynamicMutation)
                //{
                //    //==============Synamic Mutation==================
                //    for (i = 0; i < E.Parent.Parent.RowCount; i++)
                //    {
                //        Theta[i] = 0;

                //        for (p = 0; p < E.Parent.Count; p++)
                //        {
                //            Theta[i] += (E.Parent[p].RowsCoverage[i] == 1) ? (1) : (0);
                //        }

                //        if (Theta[i] < eps * E.Parent.Count)
                //        {
                //            ChangeIndex.Clear();
                //            R = new MyRandom();
                //            for (int k = 0; k < Math.Min(Md, E.Parent.Parent.Alpha[i].Count); k++)
                //            {
                //                int ind;
                //                do
                //                {
                //                    ind = E.Parent.Parent.Alpha[i][R.Next(E.Parent.Parent.Alpha[i].Count)];
                //                    //ind = 
                //                } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

                //                ChangeIndex.Add(ind);
                //                E[ind] = 1;
                //            }
                //        }
                //    }
                //}
            }

            public MemoryPopulationElement BinaryTournamentSelection()
            {
                MyRandom R = new MyRandom();

                MemoryPopulationElement t1 = _memorypopulation[R.Next(Count)];
                MemoryPopulationElement t2 = _memorypopulation[R.Next(Count)];

                MemoryPopulationElement P = (t1.fitness < t2.fitness) ? (t1) : (t2);

                return P;
            }

            public MemoryPopulationElement[] GetParents()
            {
                MemoryPopulationElement P1 = BinaryTournamentSelection();
                if (P1.unfitness == 0) return new MemoryPopulationElement[2] { P1, BinaryTournamentSelection() };
                List<MemoryPopulationElement> P2 = new List<MemoryPopulationElement>();

                List<int> RowIndexes = Parent.RowIndexes;
                //int i;

                List<int> Rp1 = P1.CoveredRows;
                List<int> Rp2;

                int Max = 0;
                foreach (MemoryPopulationElement p in _memorypopulation)
                {
                    if (p != P1)
                    {
                        Rp2 = p.CoveredRows;

                        List<int> Union = Rp1.Union(Rp2).ToList();
                        List<int> Intersection = Rp1.Intersect(Rp2).ToList();

                        int diff = Union.Count() - Intersection.Count();
                        if (diff > Max)
                        {
                            Max = diff;
                            P2.Clear();
                            P2.Add(p);
                        }
                        else
                            if (diff == Max)
                            {
                                P2.Add(p);
                            }
                    }
                }

                return new MemoryPopulationElement[2] { P1, P2.OrderBy(p => p.fitness).ToList()[0] };
            }

            public void GetGroups(MemoryPopulationElement Individum, ref List<int> GroupA, ref List<int> GroupB, ref List<int> GroupC, ref List<int> GroupD)
            {
                GroupA = new List<int>();
                GroupB = new List<int>();
                GroupC = new List<int>();
                GroupD = new List<int>();

                for (int i = 0; i < Count; i++)
                {

                    if (_memorypopulation[i].unfitness > Individum.unfitness)
                    {
                        if (_memorypopulation[i].fitness > Individum.fitness)
                        {
                            GroupA.Add(i);
                        }
                        else
                            if (_memorypopulation[i].fitness <= Individum.fitness)
                            {
                                GroupB.Add(i);
                            }
                    }
                    else
                        if (_memorypopulation[i].unfitness >= Individum.unfitness)
                        {
                            if (_memorypopulation[i].fitness > Individum.fitness)
                            {
                                GroupC.Add(i);
                            }
                            else
                                if (_memorypopulation[i].fitness <= Individum.fitness)
                                {
                                    GroupD.Add(i);
                                }
                        }
                }
            }

            public bool Exists(MemoryPopulationElement P)
            {
                return (from p in _memorypopulation where p != null && p.fitness == P.fitness && p.unfitness == P.unfitness/*p.x.SequenceEqual(P.x)*/ select p).Count() > 0;
            }

            public virtual int Replace(MemoryPopulationElement P)
            {
                List<int> GroupA = null,
                          GroupB = null,
                          GroupC = null,
                          GroupD = null,
                          Group = null;

                GetGroups(P, ref GroupA, ref GroupB, ref GroupC, ref GroupD);

                double MaxFitness;
                int MFIndex = -1;
                if (GroupA.Count > 0)
                {
                    Group = GroupA;
                    MaxFitness = (from ind in Group select _memorypopulation[ind].fitness).Max();
                    MFIndex = (from ind in Group where _memorypopulation[ind].fitness == MaxFitness orderby _memorypopulation[ind].unfitness descending select ind).ToList()[0];
                    _memorypopulation[MFIndex] = P;
                }
                else
                    if (GroupB.Count > 0)
                    {
                        //Func<int, Void> F = new Func<int,void>(
                        Group = GroupB;
                        MaxFitness = (from ind in Group select _memorypopulation[ind].fitness).Max();
                        MFIndex = (from ind in Group where _memorypopulation[ind].fitness == MaxFitness orderby _memorypopulation[ind].unfitness descending select ind).ToList()[0];
                        _memorypopulation[MFIndex] = P;
                    }
                    else
                        if (GroupC.Count > 0)
                        {
                            //if ((from ind in GroupC orderby _population[ind].unfitness descending select _population[ind].unfitness).ToList()[0] == P.unfitness)
                            {
                                Group = GroupC;
                                MaxFitness = (from ind in Group select _memorypopulation[ind].fitness).Max();
                                MFIndex = (from ind in Group where _memorypopulation[ind].fitness == MaxFitness orderby _memorypopulation[ind].unfitness descending select ind).ToList()[0];
                                _memorypopulation[MFIndex] = P;
                            }
                        }
                        else
                        {
                            //Group = GroupD;
                        }

                Group = null;
                if (Group != null)
                {
                    MaxFitness = (from ind in Group select _memorypopulation[ind].unfitness).Max();
                    MFIndex = (from ind in Group where _memorypopulation[ind].unfitness == MaxFitness orderby _memorypopulation[ind].fitness descending select ind).ToList()[0];
                    if (MaxFitness == 0)
                    {
                        //MessageBox.Show("Unfitness = 0!");
                        return -1;
                    }
                    _memorypopulation[MFIndex] = P;
                }

                return MFIndex;
            }

            public static void HeuristicFeasibilityOperator(ref MemoryPopulationElement P)
            {
                int ColCount = P.Parent.Parent.ColCount;
                int RowCount = P.Parent.Parent.RowCount;

                List<int> I = P.Parent.Parent.RowIndexes;
                List<int> J = P.Parent.Parent.ColIndexes;
                List<int>[] Alpha = P.Parent.Parent.Alpha;
                List<int>[] Beta = P.Parent.Parent.Beta;

                MemoryPopulationElement TmpP = P;
                List<int> S = (from jj in J where TmpP[jj] == 1 select jj).ToList();
                //List<int> W = P.CoveredRows;
                List<int> U, V = new List<int>();
                int[] w;// = P.RowsCoverage;
                int i, j;
                //==================Identify over-cross rows======================
                w = (from ii in I select Alpha[ii].Intersect(S).Count()).ToArray();

                List<int> T = new List<int>();
                int g = 0;
                do
                {
                    T.Clear();
                    T.AddRange(S);
                    MyRandom R = new MyRandom();
                    if (T.Count > 0)
                    {
                        do
                        {
                            j = T[R.Next(T.Count)];
                            T.Remove(j);

                            if ((from ii in Beta[j] where w[ii] >= 2 select ii).Count() == Beta[j].Count)
                            {
                                S.Remove(j);
                                foreach (int ii in Beta[j])
                                {
                                    w[ii]--;
                                }
                            }
                        } while (T.Count != 0);
                    }
                    //==================Identify under-cross rows=====================
                    U = (from ii in I where w[ii] == 0 select ii).ToList();
                    if (U.Count > 0)
                    {
                        V.Clear();
                        V.AddRange(U);

                        R = new MyRandom();
                        do
                        {
                            i = V[R.Next(V.Count)];
                            V.Remove(i);

                            var Res = (from jj in Alpha[i] where Beta[jj].Count != 0 && Beta[jj].Except(Beta[jj].Intersect(U)).Count() == 0/* U.Except(Beta[jj].Intersect(U)).Count() == 0 */ select new { jj, cc = TmpP.Parent.Parent.c[jj] / Beta[jj].Count }).OrderBy(a => a.cc);

                            if (Res.Count() > 0)
                            {
                                j = Res.ToList()[0].jj;
                                S.Add(j);
                                foreach (int ii in Beta[j])
                                {
                                    w[ii] = 1;
                                }

                                U = U.Except(Beta[j]).ToList();
                                V = V.Except(Beta[j]).ToList();
                            }
                        } while (V.Count != 0);
                    }

                    g++;
                } while (g < 1);

                for (i = 0; i < ColCount; i++)
                {
                    if (S.Contains(i))
                    {
                        P[i] = 1;
                    }
                    else
                    {
                        P[i] = 0;
                    }
                }

                P.Process();
            }

            System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
            {
                return MemoryPopulation.GetEnumerator();
            }

            public IEnumerator<MemoryPopulationElement> GetEnumerator()
            {
                foreach (MemoryPopulationElement p in _memorypopulation)
                {
                    yield return p;
                }
            }
            #endregion
        }
        #endregion
        #region ReferenceSet
        class MemoryReferenceSet : MemoryPopulationClass
        {
            #region Methods
            public MemoryReferenceSet(SPPClass Parent, int Count = 0) : base(Parent, Count) { }
            public static MemoryReferenceSet DiversificationAlgorithm(MemoryPopulationClass CL, int N)
            {
                if (N >= CL.Count) return null;

                MemoryReferenceSet RefSet = new MemoryReferenceSet(CL.Parent);

                int i;
                for (i = 0; i < N / 2; i++)
                {
                    double minUnfitness = CL.OrderBy(p => p.unfitness).ToList()[0].unfitness;
                    MemoryPopulationElement c = (from p in CL where p.unfitness == minUnfitness select p).OrderBy(p => p.fitness).ToList()[0];
                    RefSet.Add(c);
                    CL.Delete(c);
                    /*if (RefSet.Parent.UsingHeuristicOperator)
                        PopulationClass.HeuristicFeasibilityOperator(ref c);*/
                }

                for (i = 0; i < N / 2; i++)
                {
                    MemoryPopulationElement c = CL.OrderByDescending(p => p.GetDistance(RefSet)).ToList()[0];//(from p in CL select new { elem = p, dist = p.GetDistance(RefSet) }).OrderByDescending(d => d.dist).ToList()[0].elem; 
                    RefSet.Add(c);
                    CL.Delete(c);
                    /*if (RefSet.Parent.UsingHeuristicOperator)
                        PopulationClass.HeuristicFeasibilityOperator(ref c);*/
                }
                return RefSet;
            }

            public static MemoryReferenceSet PathRelinking(MemoryPopulationElement p_start, MemoryPopulationElement p_end)
            {
                MemoryReferenceSet RS = new MemoryReferenceSet(p_start.Parent.Parent);

                int i;
                MemoryPopulationElement p0 = p_start.Copy();
                double dist = p0.GetDistance(p_end);
                while (dist > 0.1)
                {
                    MemoryPopulationElement c;// = p0.Copy();
                    List<MemoryPopulationElement> Loc = new List<MemoryPopulationElement>();
                    for (i = 0; i < p_start.Count; i++)
                    {
                        c = p0.Copy();
                        c[i] = (byte)(1 - c[i]);
                        c.Process(true);
                        if (c.GetDistance(p_end) < dist)
                            Loc.Add(c);
                        //c[i] = (byte)(1 - c[i]);
                    }

                    Loc = Loc.OrderBy(p => p.Fitness).ToList();
                    i = 0;
                    while (RS.Exists(c = Loc[i]) && i < Loc.Count) i++;
                    RS.Add(c);
                    p0 = c;
                    dist = p0.GetDistance(p_end);
                }

                return RS;
            }
            #endregion
        }
        #endregion

         #region EpsPopulationElement
    public class EpsPopulationElement : System.Collections.IEnumerable
    {
        #region EpsPopulationElement Fields
        private double _Fitness;
        private double _fitness;
        private double _unfitness;
        public byte[] Value;
        private int _Count;
        private EpsPopulationClass _parent;
        private int[] _RowsCoverage;
        private List<int> _CoveredRows;
       // private MarkType _mark = MarkType.mtNew;
        //private int[] _RowsCoverage;
        #endregion

        #region EpsPopulationElement Properties
        public int[] RowsCoverage
        {
            get { return _RowsCoverage; }
            set { _RowsCoverage = value; }
        }

        public List<int> CoveredRows
        {
            get { return _CoveredRows; }
            set { _CoveredRows = value; }
        }

        public EpsPopulationClass Parent
        {
            get { return _parent; }
            set { _parent = value; }
        }

        public double Fitness
        {
            get { return _Fitness; }
            set { _Fitness = value; }
        }

        public double fitness
        {
            get { return _fitness; }
            set { _fitness = value; }
        }

        public double unfitness
        {
            get { return _unfitness; }
            set { _unfitness = value; }
        }

        public byte[] x
        {
            get { return Value; }
            set { Value = value; }
        }

        public int Count
        {
            get { return _Count; }
            set
            {
                if (value >= 0)
                {
                    _Count = value;
                    FillRandom();
                    Process();
                }
            }
        }

        public byte this[int Index]
        {
            get { return Value[Index]; }
            set { Value[Index] = value; }
        }

        public int Index
        {
            get;
            set;
        }

       // public MarkType Mark
        //{
        //    get { return _mark; }
        //    set { _mark = value; }
        //}
        #endregion

        #region EpsPopulationElement Constructors
        public EpsPopulationElement(EpsPopulationClass Parent, int Count, bool ToProcess = true)
        {
            this._parent = Parent;
            _Fitness = 0;
            _Count = Count;
            FillRandom();
            if (ToProcess) Process();
        }
        public EpsPopulationElement(EpsPopulationClass Parent, byte[]arr, int Count = 0, bool ToProcess = true)
        {
            this._parent = Parent;
            _Fitness = 0;
            _Count = Count;
            Fill(arr);
            if (ToProcess) Process();
        }

        public EpsPopulationElement(EpsPopulationClass Parent, byte[] x, bool ToProcess = true)
        {
            this._parent = Parent;
            _Fitness = 0; _Count = x.Length;
            Value = new byte[x.Length];
            x.CopyTo(Value, 0);
            if (ToProcess) Process();
        }
        public EpsPopulationElement(EpsPopulationElement p, bool FullCopy = false)
        {
            this._parent = p.Parent;
            this._Fitness = p.Fitness;
            this._unfitness = p.unfitness;
            this._fitness = p.fitness;
            this._Count = p.Count;
            //this._mark = p.Mark;
            this.Value = new byte[p.Count];

            p.x.CopyTo(this.Value, 0);

            if (FullCopy)
            {
                this._RowsCoverage = new int[Parent.Parent.RowCount];
                p.RowsCoverage.CopyTo(this._RowsCoverage, 0);

                this._CoveredRows = new List<int>();
                this._CoveredRows.AddRange(p.CoveredRows);
            }
        }
        #endregion

        #region EpsPopulationElement Methods
        public EpsPopulationElement()
        {
        }

        public EpsPopulationElement Copy(bool FullCopy = false)
        {
            return new EpsPopulationElement(this, FullCopy);
        }

        public void Process(bool GetFitnessOnly = false)
        {
            GetFitness();
            if (!GetFitnessOnly)
            {
                GetRowsCoverage();
                _CoveredRows = RowsCoverd();
            }
        }

        public void FillRandom()
        {
            MyRandom R = new MyRandom();
            Value = new byte[Count];
            for (int i = 0; i < _Count; i++)
            {
                Value[i] = (byte)R.Next(2);
            }
        }

        public void Fill(byte[] arr)
        {
         //   MyRandom R = new MyRandom();
            Value = new byte[Count];
            for (int i = 0; i < _Count; i++)
            {
                Value[i] = arr[i];
            }
        }

        private double GetFitness()
        {
            return _Fitness = SPPClass.Fitness(Value, Parent.Parent.A, Parent.Parent.c, out _fitness, out _unfitness);
        }

        public System.Collections.IEnumerator GetEnumerator()
        {
            return Value.GetEnumerator();
        }

        private void GetRowsCoverage()
        {
            _RowsCoverage = new int[Parent.Parent.RowCount];

            for (int i = 0; i < Parent.Parent.RowCount; i++)
            {
                _RowsCoverage[i] = RowCovering(i);
            }
        }

        public int RowCovering(int RowIndex)
        {
            int res = 0;
            for (int i = 0; i < Count; i++)
            {
                res += Parent.Parent.A[RowIndex][i] * this[i];
            }

            return res;
        }

        public List<int> RowsCoverd()
        {
            List<int> L = (from i in Parent.Parent.RowIndexes where RowsCoverage[i] > 0 select i).ToList();
            return L;
        }

        public double GetDistance(EpsPopulationElement Other)
        {
            if (this._Count != Other.Count)
            {
                throw new Exception("Individuals should have identical sizes!");
            }
            double Result = 0;

            for (int i = 0; i < _Count; i++)
            {
                Result += (Value[i] - Other[i]) * (Value[i] - Other[i]);
            }
            Result = Math.Sqrt(Result);

            return Result;
        }

        public double GetDistance(EpsPopulationClass Population)
        {
            return (from p in Population select this.GetDistance(p)).OrderBy(d => d).ToList()[0];
        }

        internal void SetParent(EpsPopulationClass Parent, bool ToProcess)
        {
            _parent = Parent;
            if (ToProcess)
            {
                Process(true);
            }
        }
        #endregion
    }
    #endregion
    #region EpsPopulationClass
    public class EpsPopulationClass : IEnumerable<EpsPopulationElement>
    {
        #region PopulationClass Fields
        private List<EpsPopulationElement> _epspopulation = new List<EpsPopulationElement>();

        private int _Count;
        private int _ColCount;

        private SPPClass _parent;
        #endregion

        #region PopulationClass Properties
        public SPPClass Parent
        {
            get { return _parent; }
        }

        public List<EpsPopulationElement> EpsPopulation
        {
            get { return _epspopulation; }
        }

        public int Count
        {
            get { return _Count; }
        }

        public int ColCount
        {
            get { return _ColCount; }
        }

        public EpsPopulationElement this[int Index]
        {
            get { return _epspopulation[Index]; }
            set { _epspopulation[Index] = value; }
        }

        protected EpsPopulationElement Last
        {
            get
            {
                return this[_Count - 1];
            }

            set
            {
                this[_Count - 1] = value;
            }
        }

        #endregion

        public EpsPopulationClass(SPPClass Parent, int Count = 0)
        {
            this._parent = Parent;
            _ColCount = Parent.ColCount;
            _Count = Count;
            Generate();
        }

        public void Generate(byte[] p)
        {
            MyRandom R = new MyRandom();
            //for (int k = 0; k < Ms; k++)
            //{
            //    int ind;
            //    //  do
            //    // {
            //    ind = R.Next(E.Count);
            //    //ind = 
            //    //   } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

            //    ChangeIndex.Add(ind);
            //    E[ind] = (byte)(1 - E[ind]);//(byte)R.Next(2);
            //}

            _epspopulation.Clear();// = new PopulationElement[Count];
            byte[] temp = new byte[p.Count()];
            temp = p.ToArray();
            //int k = R.Next(0, temp.Count() - 1);
            int k = (byte)R.Next(temp.Count() - 1); 
            for (int i = 0; i < Count; i++)
            {
                //k++;
                //if (k == temp.Count() -1)
                //    k--;
                //else k ++;
                //k = (byte)R.Next(temp.Count());    
                //p[k] = (byte)Math.Abs(temp[k] - 1);
                    EpsPopulationElement tmp;
                         do
                         {
                             {
                                 p[k] = (byte)Math.Abs(temp[k] - 1);
                                 tmp = new EpsPopulationElement(this, p, ColCount);
                                 if (Parent.UsingHeuristicOperator)
                                     EpsPopulationClass.HeuristicFeasibilityOperator(ref tmp);
                                 k = (byte)R.Next(temp.Count() - 1);
                                //if (k == p.Count() - 1)
                                //    k--;
                                // else k++;
                                // if (i < Count / 2)
                                //     k++;
                                // else k -= (i - 1);
                             }
                    }
                          while (this.Exists(tmp));

                    _epspopulation.Add(tmp);
                    _epspopulation[i].Index = i;
            }
        }
        public EpsPopulationClass(SPPClass Parent, int eps, byte []p)
        {
            this._parent = Parent;
            _ColCount = Parent.ColCount;
            _Count = eps;
            Generate(p);
        }

        #region EpsPopulationClass Methods
        public EpsPopulationClass Add(EpsPopulationElement p)
        {
            if (this.Exists(p)) return this;

            _epspopulation.Add(p);
            EpsPopulationElement tmp = _epspopulation[_Count];// = p.Copy();
            if (tmp.Parent.Parent == this.Parent)
                tmp.SetParent(this, false);
            else tmp.SetParent(this, true);
            _Count++;

            return this;
        }

        public EpsPopulationClass Add(EpsPopulationClass OtherPopulation)
        {
            foreach (EpsPopulationElement p in OtherPopulation)
            {
                this.Add(p);
            }

            return this;
        }

        public EpsPopulationClass Delete(EpsPopulationElement p)
        {
            _epspopulation.Remove(p);// = (from pp in _population where pp != p select pp).ToList();
            _Count--;
            return this;
        }

        public void Generate()
        {
            _epspopulation.Clear();// = new PopulationElement[Count];

            for (int i = 0; i < Count; i++)
            {
                EpsPopulationElement tmp;
                do
                {
                    tmp = new EpsPopulationElement(this, ColCount);
                    if (Parent.UsingHeuristicOperator)
                        EpsPopulationClass.HeuristicFeasibilityOperator(ref tmp);
                }
                while (this.Exists(tmp));

                _epspopulation.Add(tmp);
                _epspopulation[i].Index = i;
            }
        }

        public double Fitness(int Index)
        {
            if (Index >= 0 && Index < Count)
                return EpsPopulation[Index].Fitness;
            else return Double.MaxValue;
        }

        public void Sort()
        {
            _epspopulation = _epspopulation.OrderBy(element => element.Fitness).ToList();
        }

        public static EpsPopulationElement Crossover(EpsPopulationElement P1, EpsPopulationElement P2)
        {
            if (P1.Count == P2.Count && P1.Parent == P2.Parent)
            {
                int count = P1.Count;
                EpsPopulationElement Mask = new EpsPopulationElement(P1.Parent, count);
                EpsPopulationElement Child = new EpsPopulationElement(P1.Parent, count, false);

                for (int i = 0; i < count; i++)
                {
                    Child[i] = (Mask[i] == 0) ? (P1[i]) : (P2[i]);
                }

                Child.Process();
                return Child;
            }
            else return null;
        }


        public static void Mutation(ref EpsPopulationElement E, bool DynamicMutation = false)
        {
       //     int[] Theta = new int[E.Parent.Parent.RowCount];

//            double eps = 0.5;
            //int Md = 5, 
            int    Ms = 3;
            List<int> ChangeIndex = new List<int>();

            //============Static Mutation================
            MyRandom R = new MyRandom();
            for (int k = 0; k < Ms; k++)
            {
                int ind;
              //  do
               // {
                    ind = R.Next(E.Count);
                    //ind = 
             //   } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

                ChangeIndex.Add(ind);
                E[ind] = (byte)(1 - E[ind]);//(byte)R.Next(2);
            }

            E.Process();

            //if (DynamicMutation)
            //{
            //    //==============Synamic Mutation==================
            //    for (i = 0; i < E.Parent.Parent.RowCount; i++)
            //    {
            //        Theta[i] = 0;

            //        for (p = 0; p < E.Parent.Count; p++)
            //        {
            //            Theta[i] += (E.Parent[p].RowsCoverage[i] == 1) ? (1) : (0);
            //        }

            //        if (Theta[i] < eps * E.Parent.Count)
            //        {
            //            ChangeIndex.Clear();
            //            R = new MyRandom();
            //            for (int k = 0; k < Math.Min(Md, E.Parent.Parent.Alpha[i].Count); k++)
            //            {
            //                int ind;
            //                do
            //                {
            //                    ind = E.Parent.Parent.Alpha[i][R.Next(E.Parent.Parent.Alpha[i].Count)];
            //                    //ind = 
            //                } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

            //                ChangeIndex.Add(ind);
            //                E[ind] = 1;
            //            }
            //        }
            //    }
            //}
        }

        public EpsPopulationElement BinaryTournamentSelection()
        {
            MyRandom R = new MyRandom();

            EpsPopulationElement t1 = _epspopulation[R.Next(Count)];
            EpsPopulationElement t2 = _epspopulation[R.Next(Count)];

            EpsPopulationElement P = (t1.fitness < t2.fitness) ? (t1) : (t2);

            return P;
        }

        public EpsPopulationElement[] GetParents()
        {
            EpsPopulationElement P1 = BinaryTournamentSelection();
            if (P1.unfitness == 0) return new EpsPopulationElement[2] { P1, BinaryTournamentSelection() };
            List<EpsPopulationElement> P2 = new List<EpsPopulationElement>();

            List<int> RowIndexes = Parent.RowIndexes;
            //int i;

            List<int> Rp1 = P1.CoveredRows;
            List<int> Rp2;

            int Max = 0;
            foreach (EpsPopulationElement p in _epspopulation)
            {
                if (p != P1)
                {
                    Rp2 = p.CoveredRows;

                    List<int> Union = Rp1.Union(Rp2).ToList();
                    List<int> Intersection = Rp1.Intersect(Rp2).ToList();

                    int diff = Union.Count() - Intersection.Count();
                    if (diff > Max)
                    {
                        Max = diff;
                        P2.Clear();
                        P2.Add(p);
                    }
                    else
                        if (diff == Max)
                        {
                            P2.Add(p);
                        }
                }
            }

            return new EpsPopulationElement[2] { P1, P2.OrderBy(p => p.fitness).ToList()[0] };
        }

        public void GetGroups(EpsPopulationElement Individum, ref List<int> GroupA, ref List<int> GroupB, ref List<int> GroupC, ref List<int> GroupD)
        {
            GroupA = new List<int>();
            GroupB = new List<int>();
            GroupC = new List<int>();
            GroupD = new List<int>();

            for (int i = 0; i < Count; i++)
            {

                if (_epspopulation[i].unfitness > Individum.unfitness)
                {
                    if (_epspopulation[i].fitness > Individum.fitness)
                    {
                        GroupA.Add(i);
                    }
                    else
                        if (_epspopulation[i].fitness <= Individum.fitness)
                        {
                            GroupB.Add(i);
                        }
                }
                else
                    if (_epspopulation[i].unfitness >= Individum.unfitness)
                    {
                        if (_epspopulation[i].fitness > Individum.fitness)
                        {
                            GroupC.Add(i);
                        }
                        else
                            if (_epspopulation[i].fitness <= Individum.fitness)
                            {
                                GroupD.Add(i);
                            }
                    }
            }
        }

        public bool Exists(EpsPopulationElement P)
        {
            return (from p in _epspopulation where p != null && p.fitness == P.fitness && p.unfitness == P.unfitness/*p.x.SequenceEqual(P.x)*/ select p).Count() > 0;
        }

        public virtual int Replace(EpsPopulationElement P)
        {
            List<int> GroupA = null,
                      GroupB = null,
                      GroupC = null,
                      GroupD = null,
                      Group = null;

            GetGroups(P, ref GroupA, ref GroupB, ref GroupC, ref GroupD);

            double MaxFitness;
            int MFIndex = -1;
            if (GroupA.Count > 0)
            {
                Group = GroupA;
                MaxFitness = (from ind in Group select _epspopulation[ind].fitness).Max();
                MFIndex = (from ind in Group where _epspopulation[ind].fitness == MaxFitness orderby _epspopulation[ind].unfitness descending select ind).ToList()[0];
                _epspopulation[MFIndex] = P;
            }
            else
                if (GroupB.Count > 0)
                {
                    //Func<int, Void> F = new Func<int,void>(
                    Group = GroupB;
                    MaxFitness = (from ind in Group select _epspopulation[ind].fitness).Max();
                    MFIndex = (from ind in Group where _epspopulation[ind].fitness == MaxFitness orderby _epspopulation[ind].unfitness descending select ind).ToList()[0];
                    _epspopulation[MFIndex] = P;
                }
                else
                    if (GroupC.Count > 0)
                    {
                        //if ((from ind in GroupC orderby _population[ind].unfitness descending select _population[ind].unfitness).ToList()[0] == P.unfitness)
                        {
                            Group = GroupC;
                            MaxFitness = (from ind in Group select _epspopulation[ind].fitness).Max();
                            MFIndex = (from ind in Group where _epspopulation[ind].fitness == MaxFitness orderby _epspopulation[ind].unfitness descending select ind).ToList()[0];
                            _epspopulation[MFIndex] = P;
                        }
                    }
                    else
                    {
                        //Group = GroupD;
                    }

            Group = null;
            if (Group != null)
            {
                MaxFitness = (from ind in Group select _epspopulation[ind].unfitness).Max();
                MFIndex = (from ind in Group where _epspopulation[ind].unfitness == MaxFitness orderby _epspopulation[ind].fitness descending select ind).ToList()[0];
                if (MaxFitness == 0)
                {
                    //MessageBox.Show("Unfitness = 0!");
                    return -1;
                }
                _epspopulation[MFIndex] = P;
            }

            return MFIndex;
        }

        public static void HeuristicFeasibilityOperator(ref EpsPopulationElement P)
        {
            int ColCount = P.Parent.Parent.ColCount;
            int RowCount = P.Parent.Parent.RowCount;

            List<int> I = P.Parent.Parent.RowIndexes;
            List<int> J = P.Parent.Parent.ColIndexes;
            List<int>[] Alpha = P.Parent.Parent.Alpha;
            List<int>[] Beta = P.Parent.Parent.Beta;

            EpsPopulationElement TmpP = P;
            List<int> S = (from jj in J where TmpP[jj] == 1 select jj).ToList();
            //List<int> W = P.CoveredRows;
            List<int> U, V = new List<int>();
            int[] w;// = P.RowsCoverage;
            int i, j;
            //==================Identify over-cross rows======================
            w = (from ii in I select Alpha[ii].Intersect(S).Count()).ToArray();

            List<int> T = new List<int>();
            int g = 0;
            do
            {
                T.Clear();
                T.AddRange(S);
                MyRandom R = new MyRandom();
                if (T.Count > 0)
                {
                    do
                    {
                        j = T[R.Next(T.Count)];
                        T.Remove(j);

                        if ((from ii in Beta[j] where w[ii] >= 2 select ii).Count() == Beta[j].Count)
                        {
                            S.Remove(j);
                            foreach (int ii in Beta[j])
                            {
                                w[ii]--;
                            }
                        }
                    } while (T.Count != 0);
                }
                //==================Identify under-cross rows=====================
                U = (from ii in I where w[ii] == 0 select ii).ToList();
                if (U.Count > 0)
                {
                    V.Clear();
                    V.AddRange(U);

                    R = new MyRandom();
                    do
                    {
                        i = V[R.Next(V.Count)];
                        V.Remove(i);

                        var Res = (from jj in Alpha[i] where Beta[jj].Count != 0 && Beta[jj].Except(Beta[jj].Intersect(U)).Count() == 0/* U.Except(Beta[jj].Intersect(U)).Count() == 0 */ select new { jj, cc = TmpP.Parent.Parent.c[jj] / Beta[jj].Count }).OrderBy(a => a.cc);

                        if (Res.Count() > 0)
                        {
                            j = Res.ToList()[0].jj;
                            S.Add(j);
                            foreach (int ii in Beta[j])
                            {
                                w[ii] = 1;
                            }

                            U = U.Except(Beta[j]).ToList();
                            V = V.Except(Beta[j]).ToList();
                        }
                    } while (V.Count != 0);
                }

                g++;
            } while (g < 1);

            for (i = 0; i < ColCount; i++)
            {
                if (S.Contains(i))
                {
                    P[i] = 1;
                }
                else
                {
                    P[i] = 0;
                }
            }

            P.Process();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return EpsPopulation.GetEnumerator();
        }

        public IEnumerator<EpsPopulationElement> GetEnumerator()
        {
            foreach (EpsPopulationElement p in _epspopulation)
            {
                yield return p;
            }
        }
        #endregion
    }
    #endregion
    #region ReferenceSet
    class EpsReferenceSet : EpsPopulationClass
    {
        #region Methods
        public EpsReferenceSet(SPPClass Parent, int Count = 0) : base(Parent, Count) { }
        public static EpsReferenceSet DiversificationAlgorithm(EpsPopulationClass CL, int N)
        {
            if (N >= CL.Count) return null;

            EpsReferenceSet RefSet = new EpsReferenceSet(CL.Parent);

            int i;
            for (i = 0; i < N / 2; i++)
            {
                double minUnfitness = CL.OrderBy(p => p.unfitness).ToList()[0].unfitness;
                EpsPopulationElement c = (from p in CL where p.unfitness == minUnfitness select p).OrderBy(p => p.fitness).ToList()[0];
                RefSet.Add(c);
                CL.Delete(c);
                /*if (RefSet.Parent.UsingHeuristicOperator)
                    PopulationClass.HeuristicFeasibilityOperator(ref c);*/
            }

            for (i = 0; i < N / 2; i++)
            {
                EpsPopulationElement c = CL.OrderByDescending(p => p.GetDistance(RefSet)).ToList()[0];//(from p in CL select new { elem = p, dist = p.GetDistance(RefSet) }).OrderByDescending(d => d.dist).ToList()[0].elem; 
                RefSet.Add(c);
                CL.Delete(c);
                /*if (RefSet.Parent.UsingHeuristicOperator)
                    PopulationClass.HeuristicFeasibilityOperator(ref c);*/
            }
            return RefSet;
        }

        public static EpsReferenceSet PathRelinking(EpsPopulationElement p_start, EpsPopulationElement p_end)
        {
            EpsReferenceSet RS = new EpsReferenceSet(p_start.Parent.Parent);

            int i;
            EpsPopulationElement p0 = p_start.Copy();
            double dist = p0.GetDistance(p_end);
            while (dist > 0.1)
            {
                EpsPopulationElement c;// = p0.Copy();
                List<EpsPopulationElement> Loc = new List<EpsPopulationElement>();
                for (i = 0; i < p_start.Count; i++)
                {
                    c = p0.Copy();
                    c[i] = (byte)(1 - c[i]);
                    c.Process(true);
                    if (c.GetDistance(p_end) < dist)
                        Loc.Add(c);
                    //c[i] = (byte)(1 - c[i]);
                }

                Loc = Loc.OrderBy(p => p.Fitness).ToList();
                i = 0;
                while (RS.Exists(c = Loc[i]) && i < Loc.Count) i++;
                RS.Add(c);
                p0 = c;
                dist = p0.GetDistance(p_end);
            }

            return RS;
        }
        #endregion
    }
    #endregion

        #region PopulationElement
        public class PopulationElement : System.Collections.IEnumerable
        {
            #region PopulationElement Fields
            private double _Fitness;
            private double _fitness;
            private double _unfitness;
            public byte[] Value;
            private int _Count;
            private PopulationClass _parent;
            private int[] _RowsCoverage;
            private List<int> _CoveredRows;
            private MarkType _mark = MarkType.mtNew;
            //private int[] _RowsCoverage;
            #endregion

            #region PopulationElement Properties
            public int[] RowsCoverage
            {
                get { return _RowsCoverage; }
                set { _RowsCoverage = value; }
            }

            public List<int> CoveredRows
            {
                get { return _CoveredRows; }
                set { _CoveredRows = value; }
            }

            public PopulationClass Parent
            {
                get { return _parent; }
            }

            public double Fitness
            {
                get { return _Fitness; }
                set { _Fitness = value; }
            }

            public double fitness
            {
                get { return _fitness; }
                set { _fitness = value; }
            }

            public double unfitness
            {
                get { return _unfitness; }
                set { _unfitness = value; }
            }

            public byte[] x
            {
                get { return Value; }
                set { Value = value; }
            }

            public int Count
            {
                get { return _Count; }
                set
                {
                    if (value >= 0)
                    {
                        _Count = value;
                        FillRandom();
                        Process();
                    }
                }
            }

            public byte this[int Index]
            {
                get { return Value[Index]; }
                set { Value[Index] = value; }
            }

            public int Index
            {
                get;
                set;
            }

            public MarkType Mark
            {
                get { return _mark; }
                set { _mark = value; }
            }
            #endregion

            #region PopulationElement Constructors
            public PopulationElement(PopulationClass Parent, int Count, bool ToProcess = true)
            {
                this._parent = Parent;
                _Fitness = 0;
                _Count = Count;
                FillRandom();
                if (ToProcess) Process();
            }

            public PopulationElement(PopulationClass Parent, byte[] x, bool ToProcess = true)
            {
                this._parent = Parent;
                _Fitness = 0; _Count = x.Length;
                Value = new byte[x.Length];
                x.CopyTo(Value, 0);
                if (ToProcess) Process();
            }
            public PopulationElement(PopulationElement p, bool FullCopy = false)
            {
                this._parent = p.Parent;
                this._Fitness = p.Fitness;
                this._unfitness = p.unfitness;
                this._fitness = p.fitness;
                this._Count = p.Count;
                this._mark = p.Mark;
                this.Value = new byte[p.Count];

                p.x.CopyTo(this.Value, 0);

                if (FullCopy)
                {
                    this._RowsCoverage = new int[Parent.Parent.RowCount];
                    p.RowsCoverage.CopyTo(this._RowsCoverage, 0);

                    this._CoveredRows = new List<int>();
                    this._CoveredRows.AddRange(p.CoveredRows);
                }
            }
            #endregion

            #region PopulationElement Methods
            public PopulationElement() { }

            public PopulationElement Convert(EpsPopulationElement eps, PopulationElement p)
            {
                p.x = eps.x;
                p.unfitness = eps.unfitness;
                p.RowsCoverage = eps.RowsCoverage;
                p.Fitness = eps.Fitness;
                p.fitness = eps.fitness;
                p.CoveredRows = eps.CoveredRows;
               // p.Count = eps.x.Count();
               // p.unfitness = eps.unfitness;
                return p;
            }

            public PopulationElement Convert(MemoryPopulationElement memory, PopulationElement p)
            {
                memory.x = p.x;
                memory.unfitness = p.unfitness;
                memory.RowsCoverage = p.RowsCoverage;
                memory.Fitness = p.Fitness;
                memory.fitness = p.fitness;
                memory.CoveredRows = p.CoveredRows;
                // p.Count = eps.x.Count();
                // p.unfitness = eps.unfitness;
                return p;
            }

            public PopulationElement Copy(bool FullCopy = false)
            {
                return new PopulationElement(this, FullCopy);
            }

            public void Process(bool GetFitnessOnly = false)
            {
                GetFitness();
                if (!GetFitnessOnly)
                {
                    GetRowsCoverage();
                    _CoveredRows = RowsCoverd();
                }
            }

            public void FillRandom()
            {
                MyRandom R = new MyRandom();
                Value = new byte[Count];
                for (int i = 0; i < _Count; i++)
                {
                    Value[i] = (byte)R.Next(2);
                }
            }

            private double GetFitness()
            {
                return _Fitness = SPPClass.Fitness(Value, Parent.Parent.A, Parent.Parent.c, out _fitness, out _unfitness);
            }

            public System.Collections.IEnumerator GetEnumerator()
            {
                return Value.GetEnumerator();
            }

            private void GetRowsCoverage()
            {
                _RowsCoverage = new int[Parent.Parent.RowCount];

                for (int i = 0; i < Parent.Parent.RowCount; i++)
                {
                    _RowsCoverage[i] = RowCovering(i);
                }
            }

            public int RowCovering(int RowIndex)
            {
                int res = 0;
                for (int i = 0; i < Count; i++)
                {
                    res += Parent.Parent.A[RowIndex][i] * this[i];
                }

                return res;
            }

            public List<int> RowsCoverd()
            {
                List<int> L = (from i in Parent.Parent.RowIndexes where RowsCoverage[i] > 0 select i).ToList();
                return L;
            }

            public double GetDistance(PopulationElement Other)
            {
                if (this._Count != Other.Count)
                {
                    throw new Exception("Individuals should have identical sizes!");
                }
                double Result = 0;

                for (int i = 0; i < _Count; i++)
                {
                    Result += (Value[i] - Other[i]) * (Value[i] - Other[i]);
                }
                Result = Math.Sqrt(Result);

                return Result;
            }

            public double GetDistance(PopulationClass Population)
            {
                return (from p in Population select this.GetDistance(p)).OrderBy(d => d).ToList()[0];
            }

            internal void SetParent(PopulationClass Parent, bool ToProcess)
            {
                _parent = Parent;
                if (ToProcess)
                {
                    Process(true);
                }
            }
            #endregion
        }
        #endregion
        #region PopulationClass
        public class PopulationClass : IEnumerable<PopulationElement>
        {
            #region PopulationClass Fields
            private List<PopulationElement> _population = new List<PopulationElement>();

            private int _Count;
            private int _ColCount;

            private SPPClass _parent;
            #endregion

            #region PopulationClass Properties
            public SPPClass Parent
            {
                get { return _parent; }
            }

            public List<PopulationElement> Population
            {
                get { return _population; }
            }

            public int Count
            {
                get { return _Count; }
            }

            public int ColCount
            {
                get { return _ColCount; }
            }

            public PopulationElement this[int Index]
            {
                get { return _population[Index]; }
                set { _population[Index] = value; }
            }

            protected PopulationElement Last
            {
                get
                {
                    return this[_Count - 1];
                }

                set
                {
                    this[_Count - 1] = value;
                }
            }

            #endregion

            public PopulationClass(SPPClass Parent, int Count, PopulationElement p)
            {
                this._parent = Parent;
                _ColCount = Parent.ColCount;
                _Count = Count;
                Generate(p);
            }

            public PopulationClass(SPPClass Parent, int Count = 0)
            {
                this._parent = Parent;
                _ColCount = Parent.ColCount;
                _Count = Count;
                Generate();
            }

            #region PopulationClass Methods
            public PopulationClass Add(PopulationElement p)
            {
                if (this.Exists(p)) return this;

                _population.Add(p);
                PopulationElement tmp = _population[_Count];// = p.Copy();
                if (tmp.Parent.Parent == this.Parent)
                    tmp.SetParent(this, false);
                else tmp.SetParent(this, true);
                _Count++;

                return this;
            }

            public PopulationClass Add(PopulationClass OtherPopulation)
            {
                foreach (PopulationElement p in OtherPopulation)
                {
                    this.Add(p);
                }

                return this;
            }

            public PopulationClass Delete(PopulationElement p)
            {
                _population.Remove(p);// = (from pp in _population where pp != p select pp).ToList();
                _Count--;
                return this;
            }

            public void Generate()
            {
                _population.Clear();// = new PopulationElement[Count];

                for (int i = 0; i < Count; i++)
                {
                    PopulationElement tmp;
                    do
                    {
                        tmp = new PopulationElement(this, ColCount);
                        if (Parent.UsingHeuristicOperator)
                            PopulationClass.HeuristicFeasibilityOperator(ref tmp);
                    }
                    while (this.Exists(tmp));

                    _population.Add(tmp);
                    _population[i].Index = i;
                }
            }
            public void Generate(PopulationElement tmp)
            {
                _population.Clear();// = new PopulationElement[Count];

                for (int i = 0; i < Count; i++)
                {
                    do
                    {
                        // tmp = new PopulationElement(this, ColCount);
                        PopulationClass.Mutation(ref tmp);
                        if (Parent.UsingHeuristicOperator)
                            PopulationClass.HeuristicFeasibilityOperator(ref tmp);
                    }
                    while (this.Exists(tmp));

                    _population.Add(tmp);
                    _population[i].Index = i;
                }
            }

            public double Fitness(int Index)
            {
                if (Index >= 0 && Index < Count)
                    return Population[Index].Fitness;
                else return Double.MaxValue;
            }

            public void Sort()
            {
                _population = _population.OrderBy(element => element.Fitness).ToList();
            }

            public static PopulationElement Crossover(PopulationElement P1, PopulationElement P2)
            {
                if (P1.Count == P2.Count && P1.Parent == P2.Parent)
                {
                    int count = P1.Count;
                    PopulationElement Mask = new PopulationElement(P1.Parent, count);
                    PopulationElement Child = new PopulationElement(P1.Parent, count, false);

                    for (int i = 0; i < count; i++)
                    {
                        Child[i] = (Mask[i] == 0) ? (P1[i]) : (P2[i]);
                    }

                    Child.Process();
                    return Child;
                }
                else return null;
            }


            public static void Mutation(ref PopulationElement E, bool DynamicMutation = false)
            {
                int[] Theta = new int[E.Parent.Parent.RowCount];

                int i, p;
                double eps = 0.5;
                int Md = 5, Ms = 3;
                List<int> ChangeIndex = new List<int>();

                //============Static Mutation================
                MyRandom R = new MyRandom();
                for (int k = 0; k < Ms; k++)
                {
                    int ind;
                    do
                    {
                        ind = R.Next(E.Count);
                        //ind = 
                    } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

                    ChangeIndex.Add(ind);
                    E[ind] = (byte)(1 - E[ind]);//(byte)R.Next(2);
                }

                E.Process();

                if (DynamicMutation)
                {
                    //==============Synamic Mutation==================
                    for (i = 0; i < E.Parent.Parent.RowCount; i++)
                    {
                        Theta[i] = 0;

                        for (p = 0; p < E.Parent.Count; p++)
                        {
                            Theta[i] += (E.Parent[p].RowsCoverage[i] == 1) ? (1) : (0);
                        }

                        if (Theta[i] < eps * E.Parent.Count)
                        {
                            ChangeIndex.Clear();
                            R = new MyRandom();
                            for (int k = 0; k < Math.Min(Md, E.Parent.Parent.Alpha[i].Count); k++)
                            {
                                int ind;
                                do
                                {
                                    ind = E.Parent.Parent.Alpha[i][R.Next(E.Parent.Parent.Alpha[i].Count)];
                                    //ind = 
                                } while ((from val in ChangeIndex where val == ind select val).Count() != 0);

                                ChangeIndex.Add(ind);
                                E[ind] = 1;
                            }
                        }
                    }
                }

                E.Process();
            }

            public PopulationElement BinaryTournamentSelection()
            {
                MyRandom R = new MyRandom();

                PopulationElement t1 = _population[R.Next(Count)];
                PopulationElement t2 = _population[R.Next(Count)];

                PopulationElement P = (t1.fitness < t2.fitness) ? (t1) : (t2);

                return P;
            }

            public PopulationElement[] GetParents()
            {
                PopulationElement P1 = BinaryTournamentSelection();
                if (P1.unfitness == 0) return new PopulationElement[2] { P1, BinaryTournamentSelection() };
                List<PopulationElement> P2 = new List<PopulationElement>();

                List<int> RowIndexes = Parent.RowIndexes;
                //int i;

                List<int> Rp1 = P1.CoveredRows;
                List<int> Rp2;

                int Max = 0;
                foreach (PopulationElement p in _population)
                {
                    if (p != P1)
                    {
                        Rp2 = p.CoveredRows;

                        List<int> Union = Rp1.Union(Rp2).ToList();
                        List<int> Intersection = Rp1.Intersect(Rp2).ToList();

                        int diff = Union.Count() - Intersection.Count();
                        if (diff > Max)
                        {
                            Max = diff;
                            P2.Clear();
                            P2.Add(p);
                        }
                        else
                            if (diff == Max)
                            {
                                P2.Add(p);
                            }
                    }
                }

                return new PopulationElement[2] { P1, P2.OrderBy(p => p.fitness).ToList()[0] };
            }

            public void GetGroups(PopulationElement Individum, ref List<int> GroupA, ref List<int> GroupB, ref List<int> GroupC, ref List<int> GroupD)
            {
                GroupA = new List<int>();
                GroupB = new List<int>();
                GroupC = new List<int>();
                GroupD = new List<int>();

                for (int i = 0; i < Count; i++)
                {

                    if (_population[i].unfitness > Individum.unfitness)
                    {
                        if (_population[i].fitness > Individum.fitness)
                        {
                            GroupA.Add(i);
                        }
                        else
                            if (_population[i].fitness <= Individum.fitness)
                            {
                                GroupB.Add(i);
                            }
                    }
                    else
                        if (_population[i].unfitness >= Individum.unfitness)
                        {
                            if (_population[i].fitness > Individum.fitness)
                            {
                                GroupC.Add(i);
                            }
                            else
                                if (_population[i].fitness <= Individum.fitness)
                                {
                                    GroupD.Add(i);
                                }
                        }
                }
            }

            public bool Exists(PopulationElement P)
            {
                return (from p in _population where p != null && p.fitness == P.fitness && p.unfitness == P.unfitness/*p.x.SequenceEqual(P.x)*/ select p).Count() > 0;
            }

            public virtual int Replace(PopulationElement P)
            {
                List<int> GroupA = null,
                          GroupB = null,
                          GroupC = null,
                          GroupD = null,
                          Group = null;

                GetGroups(P, ref GroupA, ref GroupB, ref GroupC, ref GroupD);

                double MaxFitness;
                int MFIndex = -1;
                if (GroupA.Count > 0)
                {
                    Group = GroupA;
                    MaxFitness = (from ind in Group select _population[ind].fitness).Max();
                    MFIndex = (from ind in Group where _population[ind].fitness == MaxFitness orderby _population[ind].unfitness descending select ind).ToList()[0];
                    _population[MFIndex] = P;
                }
                else
                    if (GroupB.Count > 0)
                    {
                        //Func<int, Void> F = new Func<int,void>(
                        Group = GroupB;
                        MaxFitness = (from ind in Group select _population[ind].fitness).Max();
                        MFIndex = (from ind in Group where _population[ind].fitness == MaxFitness orderby _population[ind].unfitness descending select ind).ToList()[0];
                        _population[MFIndex] = P;
                    }
                    else
                        if (GroupC.Count > 0)
                        {
                            //if ((from ind in GroupC orderby _population[ind].unfitness descending select _population[ind].unfitness).ToList()[0] == P.unfitness)
                            {
                                Group = GroupC;
                                MaxFitness = (from ind in Group select _population[ind].fitness).Max();
                                MFIndex = (from ind in Group where _population[ind].fitness == MaxFitness orderby _population[ind].unfitness descending select ind).ToList()[0];
                                _population[MFIndex] = P;
                            }
                        }
                        else
                        {
                            //Group = GroupD;
                        }

                Group = null;
                if (Group != null)
                {
                    MaxFitness = (from ind in Group select _population[ind].unfitness).Max();
                    MFIndex = (from ind in Group where _population[ind].unfitness == MaxFitness orderby _population[ind].fitness descending select ind).ToList()[0];
                    if (MaxFitness == 0)
                    {
                        //MessageBox.Show("Unfitness = 0!");
                        return -1;
                    }
                    _population[MFIndex] = P;
                }

                return MFIndex;
            }

            public static void HeuristicFeasibilityOperator(ref PopulationElement P)
            {
                int ColCount = P.Parent.Parent.ColCount;
                int RowCount = P.Parent.Parent.RowCount;

                List<int> I = P.Parent.Parent.RowIndexes;
                List<int> J = P.Parent.Parent.ColIndexes;
                List<int>[] Alpha = P.Parent.Parent.Alpha;
                List<int>[] Beta = P.Parent.Parent.Beta;

                PopulationElement TmpP = P;
                List<int> S = (from jj in J where TmpP[jj] == 1 select jj).ToList();
                //List<int> W = P.CoveredRows;
                List<int> U, V = new List<int>();
                int[] w;// = P.RowsCoverage;
                int i, j;
                //==================Identify over-cross rows======================
                w = (from ii in I select Alpha[ii].Intersect(S).Count()).ToArray();

                List<int> T = new List<int>();
                int g = 0;
                do
                {
                    T.Clear();
                    T.AddRange(S);
                    MyRandom R = new MyRandom();
                    if (T.Count > 0)
                    {
                        do
                        {
                            j = T[R.Next(T.Count)];
                            T.Remove(j);

                            if ((from ii in Beta[j] where w[ii] >= 2 select ii).Count() == Beta[j].Count)
                            {
                                S.Remove(j);
                                foreach (int ii in Beta[j])
                                {
                                    w[ii]--;
                                }
                            }
                        } while (T.Count != 0);
                    }
                    //==================Identify under-cross rows=====================
                    U = (from ii in I where w[ii] == 0 select ii).ToList();
                    if (U.Count > 0)
                    {
                        V.Clear();
                        V.AddRange(U);

                        R = new MyRandom();
                        do
                        {
                            i = V[R.Next(V.Count)];
                            V.Remove(i);

                            var Res = (from jj in Alpha[i] where Beta[jj].Count != 0 && Beta[jj].Except(Beta[jj].Intersect(U)).Count() == 0/* U.Except(Beta[jj].Intersect(U)).Count() == 0 */ select new { jj, cc = TmpP.Parent.Parent.c[jj] / Beta[jj].Count }).OrderBy(a => a.cc);

                            if (Res.Count() > 0)
                            {
                                j = Res.ToList()[0].jj;
                                S.Add(j);
                                foreach (int ii in Beta[j])
                                {
                                    w[ii] = 1;
                                }

                                U = U.Except(Beta[j]).ToList();
                                V = V.Except(Beta[j]).ToList();
                            }
                        } while (V.Count != 0);
                    }

                    g++;
                } while (g < 1);

                for (i = 0; i < ColCount; i++)
                {
                    if (S.Contains(i))
                    {
                        P[i] = 1;
                    }
                    else
                    {
                        P[i] = 0;
                    }
                }

                P.Process();
            }

            System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
            {
                return Population.GetEnumerator();
            }

            public IEnumerator<PopulationElement> GetEnumerator()
            {
                foreach (PopulationElement p in _population)
                {
                    yield return p;
                }
            }
            #endregion
        }
        #endregion
        #region ReferenceSet
        class ReferenceSet : PopulationClass
        {
            #region Methods
            public ReferenceSet(SPPClass Parent, int Count = 0) : base(Parent, Count) { }
            public static ReferenceSet DiversificationAlgorithm(PopulationClass CL, int N)
            {
                if (N >= CL.Count) return null;

                ReferenceSet RefSet = new ReferenceSet(CL.Parent);

                int i;
                for (i = 0; i < N / 2; i++)
                {
                    double minUnfitness = CL.OrderBy(p => p.unfitness).ToList()[0].unfitness;
                    PopulationElement c = (from p in CL where p.unfitness == minUnfitness select p).OrderBy(p => p.fitness).ToList()[0];
                    RefSet.Add(c);
                    CL.Delete(c);
                    //if (RefSet.Parent.UsingHeuristicOperator)
                    //    PopulationClass.HeuristicFeasibilityOperator(ref c);
                }

                for (i = 0; i < N / 2; i++)
                {
                    PopulationElement c = CL.OrderByDescending(p => p.GetDistance(RefSet)).ToList()[0];//(from p in CL select new { elem = p, dist = p.GetDistance(RefSet) }).OrderByDescending(d => d.dist).ToList()[0].elem; 
                    RefSet.Add(c);
                    CL.Delete(c);
                    //if (RefSet.Parent.UsingHeuristicOperator)
                    //    PopulationClass.HeuristicFeasibilityOperator(ref c);
                }
                return RefSet;
            }

            public static ReferenceSet PathRelinking(PopulationElement p_start, PopulationElement p_end)
            {
                ReferenceSet RS = new ReferenceSet(p_start.Parent.Parent);

                int i;
                PopulationElement p0 = p_start.Copy();
                double dist = p0.GetDistance(p_end);
                while (dist > 0.1)
                {
                    PopulationElement c;// = p0.Copy();
                    List<PopulationElement> Loc = new List<PopulationElement>();
                    for (i = 0; i < p_start.Count; i++)
                    {
                        c = p0.Copy();
                        c[i] = (byte)(1 - c[i]);
                        c.Process(true);
                        if (c.GetDistance(p_end) < dist)
                            Loc.Add(c);
                    }

                    Loc = Loc.OrderBy(p => p.Fitness).ToList();
                    i = 0;
                    while (RS.Exists(c = Loc[i]) && i < Loc.Count) i++;
                    RS.Add(c);
                    p0 = c;
                    dist = p0.GetDistance(p_end);
                }

                return RS;
            }
            #endregion
        }
        #endregion
        #endregion
        public class Memory
        {
            public Memory() { }
            double _fitness;
            double _unfitness;
            List<byte> Value;

            public double fitness
            {
                get { return _fitness; }
                set { _fitness = value; }
            }

            public double unfitness
            {
                get { return _unfitness; }
                set { _unfitness = value; }
            }

            public List<byte> x
            {
                get { return Value; }
                set { Value = value; }
            }
        }

        #region Fields
        public int M = 10000;
        private byte[][] _A;
        private int[] _c;
        private byte[] _x;
        private string _FileName;
        private int _RowCount;
        private int _ColCount;
        private double _min = double.MaxValue;
        private List<int> _RowIndexes = new List<int>(),
                          _ColIndexes = new List<int>();
        private List<int>[] _Alpha, _Beta;

        public PopulationClass Population;
      //  public PopulationClass X;
        public EpsPopulationClass EpsPopulation;
        public MemoryPopulationClass MemoryPopulation;
        public Chart Ch, ChartFitness, ChartUnfitness;

        public DataGridView data, data_vectors;
        public int[] sol;

        public TextBox tbIteration;
        public bool Solving = false;
        public bool DynamicMutation = false;
        public TextBox tbOptimum;
        public bool SleepTime = false;
        public bool UsingHeuristicOperator = true;
        #endregion

        #region Properties
        public byte[][] A
        {
            get { return _A; }
        }

        public int[] c
        {
            get { return _c; }
        }

        public byte[] x
        {
            get { return _x; }
        }

        public string FileName
        {
            get { return _FileName; }
            set
            {
                if (value != "")
                {
                    _FileName = value;
                    Open();
                }
            }
        }

        public int RowCount
        {
            get { return _RowCount; }
        }

        public int ColCount
        {
            get { return _ColCount; }
        }

        public AlgType Method
        {
            get;
            set;
        }

        public double OptimalValue
        {
            get { return _min; }
        }

        public List<int> RowIndexes
        {
            get { return _RowIndexes; }
        }

        public List<int> ColIndexes
        {
            get { return _ColIndexes; }
        }

        public List<int>[] Alpha
        {
            get { return _Alpha; }
        }

        public List<int>[] Beta
        {
            get { return _Beta; }
        }

        #endregion

        #region Constructors
        public SPPClass() : this("") { }

        public SPPClass(string FileName) : this(FileName, AlgType.atGeneticAlgorithm) { }

        public SPPClass(string FileName, AlgType Method)
        {
            _ColCount = 0;
            _RowCount = 0;
            this.Method = Method;
            this.FileName = FileName;
        }
        #endregion

        #region Methods
        public void Open()
        {
            if (_FileName == "" || !File.Exists(_FileName))
            {
                throw new FileNotFoundException("File not found!", _FileName);
            }
            else
            {
                try
                {
                    using (StreamReader sr = File.OpenText(_FileName))
                    {

                        string S;
                        while ((S = sr.ReadLine()) == "") ;

                        int[] Tokens = (from str in S.Split(' ', ';') where str != "" select int.Parse(str)).ToArray<int>();

                        _RowCount = Tokens[0];
                        _ColCount = Tokens[1];

                        int i, j;
                        _RowIndexes.Clear();
                        _ColIndexes.Clear();
                        _Alpha = new List<int>[_RowCount];
                        _Beta = new List<int>[_ColCount];

                        for (i = 0; i < _RowCount; i++)
                        {
                            _RowIndexes.Add(i);
                            _Alpha[i] = new List<int>();
                        }

                        for (j = 0; j < _ColCount; j++)
                        {
                            _ColIndexes.Add(j);
                            _Beta[j] = new List<int>();
                        }

                        _A = new byte[_RowCount][];
                        for (i = 0; i < _RowCount; i++)
                        {
                            _A[i] = new byte[_ColCount];
                        }
                        _x = new byte[_ColCount];
                        _c = new int[_ColCount];

                        i = -1;
                        while ((S = sr.ReadLine()) != null)
                        {
                            if (S.Trim() == "") continue;
                            i++;
                            Tokens = (from str in S.Split(' ', ';') where str != "" select int.Parse(str)).ToArray<int>();

                            _c[i] = Tokens[0];
                            int k = Tokens[1];
                            for (j = 0; j < k; j++)
                            {
                                _A[Tokens[j + 2] - 1][i] = 1;
                                _Beta[i].Add(Tokens[j + 2] - 1);
                                _Alpha[Tokens[j + 2] - 1].Add(i);
                            }
                        }
                    }

                    _min = double.MaxValue;
                }
                catch (IOException)
                {
                    System.Windows.Forms.MessageBox.Show(String.Format("An error occured while reading file \"{0}\"", _FileName), "SPPClass: I/O Error");
                }
            }
        }

        public void Open(string FileName)
        {
            this.FileName = FileName;
        }

        public void Solve()
        {
            Solving = true; int eps = 45;
            DateTime T0 = DateTime.Now;
            int dataColumnCount = data_vectors.ColumnCount;
            switch (Method)
            {
                case AlgType.atScatterSearch:
                    {
                        int RefSetSize = 25;
                        int InitialCount =100;
                        PopulationClass CL = new PopulationClass(this, InitialCount);
                        ReferenceSet RefSet = ReferenceSet.DiversificationAlgorithm(CL, RefSetSize);

                        #region Drawing
                        if (Ch != null)
                        {
                            Ch.Parent.Invoke(new VoidFunc(Ch.Series[0].Points.Clear), new object[] { });
                            foreach (PopulationElement p in RefSet)
                            {
                                Ch.Parent.Invoke(new Func<double, double, int>(Ch.Series[0].Points.AddXY), new object[] { p.unfitness, p.fitness });
                            }
                        }

                        if (ChartUnfitness != null)
                        {
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[0].Points.Clear), new object[] { });
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[1].Points.Clear), new object[] { });
                        }

                        if (ChartFitness != null)
                        {
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[0].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[1].Points.Clear), new object[] { });
                        }
                        for (int j = 0; j < dataColumnCount; j++)
                        {
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Remove(Convert.ToString(j)); }));
                            data.Invoke(new Action(delegate { data.Columns.Remove(Convert.ToString(j)); }));
                        }
                        data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                        data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                        data_vectors.Invoke(new Action(delegate { data_vectors.Rows.Add(ColCount - 1); }));
                        data.Invoke(new Action(delegate { data.Rows.Add(RowCount - 1); }));
                        #endregion
                        int i = 1;
                        do
                        {
                            ReferenceSet Pool = new ReferenceSet(this);
                            List<PopulationElement> NewIndividuals = (from p in RefSet where p.Mark == MarkType.mtNew select p).ToList();
                            bool doPathRelinking = !(NewIndividuals.Count == RefSet.Count);
                            foreach (PopulationElement p1 in NewIndividuals)
                            {
                                foreach (PopulationElement p2 in RefSet)
                                {
                                    if (p1 != p2)
                                    {
                                        PopulationElement Child = PopulationClass.Crossover(p1, p2);
                                        PopulationClass.Mutation(ref Child, DynamicMutation);
                                        if (UsingHeuristicOperator) PopulationClass.HeuristicFeasibilityOperator(ref Child);
                                        Pool.Add(Child);
                                        if (doPathRelinking && Pool.Count < InitialCount * 2) //&& (new MyRandom()).Next(3) == 1)
                                            Pool.Add(ReferenceSet.PathRelinking(p1, p2));
                                        //Pool[Pool.Count-1].Mark=MarkType.mtNew;
                                    }
                                }
                            }

                            foreach (PopulationElement p in RefSet) p.Mark = MarkType.mtOld;

                            CL = new PopulationClass(this);
                            CL.Add(RefSet).Add(Pool);
                            RefSet = ReferenceSet.DiversificationAlgorithm(CL, RefSetSize);
                            if (RefSet == null) break;
                            #region Drawing
                            if (Ch != null)
                            {
                                Ch.Parent.Invoke(new VoidFunc(Ch.Series[0].Points.Clear), new object[] { });
                                foreach (PopulationElement p in RefSet)
                                {
                                    Ch.Parent.Invoke(new Func<double, double, int>(Ch.Series[0].Points.AddXY), new object[] { p.unfitness, p.fitness });
                                }
                            }

                            if (tbIteration != null)
                                tbIteration.Invoke(new VoidFunc<string>(text => tbIteration.Text = text), new object[] { String.Format("{0}/{1} Time: {2}", i, M, DateTime.Now.Subtract(T0)) });

                            if (ChartUnfitness != null)
                            {
                                List<PopulationElement> tmp = RefSet.Population.OrderBy(p => p.unfitness).ToList();
                                ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[0].Points.AddXY), new object[] { i - 1, new object[] { tmp[0].unfitness } });
                                ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[1].Points.AddXY), new object[] { i - 1, new object[] { tmp[tmp.Count - 1].unfitness } });
                            }

                            if (ChartFitness != null)
                            {
                                List<PopulationElement> tmp = RefSet.Population.OrderBy(p => p.fitness).ToList();
                                ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[0].Points.AddXY), new object[] { i - 1, new object[] { tmp[0].fitness } });
                                ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[1].Points.AddXY), new object[] { i - 1, new object[] { tmp[tmp.Count - 1].fitness } });
                            }
                            sol = new int[RowCount];
                            double minUnfitness = RefSet.OrderBy(p => p.unfitness).ToList()[0].unfitness;
                            PopulationElement _tmp = (from p in RefSet where p.unfitness == minUnfitness select p).OrderBy(p => p.fitness).ToList()[0];
                           // List<PopulationElement> _tmp = RefSet.Population.OrderBy(p => p.unfitness).ToList();
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                            data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                            for (int m = 0; m < ColCount; m++)
                                data_vectors[i-1, m].Value = _tmp.Value[m];
                            sol = new int[RowCount];
                            for (int n = 0; n < RowCount; n++)
                            {
                                for (int m = 0; m < ColCount; m++)
                                    sol[n] += _tmp.Value[m] * A[n][m];
                                data[i-1, n].Value = sol[n];
                            }
                            #endregion

                            PopulationElement[] PP = (from p in RefSet.Population where p.unfitness == 0 select p).ToArray();
                            if (PP.Length > 0)
                            {

                                PopulationElement Solvation = PP.OrderBy(p => p.fitness).ToList()[0];

                                if (Solvation.fitness < _min)
                                    _min = Solvation.fitness;

                                if (tbOptimum != null)
                                    tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = _min.ToString()), new object[] { "" });

                            }
                            else
                                if (_min == Double.MaxValue)
                                {
                                    if (tbOptimum != null)
                                        tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = "No Solution found"), new object[] { "" });
                                }
                            i++;
                        } while (i < M && Solving /* && ((from p in RefSet where p.Mark == MarkType.mtOld select p).Count() < RefSet.Count)/* && (from p in Population where p.unfitness==0 select p.unfitness).Count()!=Population.Count*/);
                        if (tbOptimum != null)
                        {
                            string result = tbOptimum.Text;
                            tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = "wszystko: optimum = " + result), new object[] { "" });
                        }
                        break;
                    }

                case AlgType.atTabuSearch:
                    {
                        PopulationClass X = new PopulationClass(this, 1);
                        PopulationElement xx = new PopulationElement();
                        List<PopulationElement> MemoryList = new List<PopulationElement>();

                        MemoryList.Add(X[0]);

                        #region Drawing
                        //if (Ch != null)
                        //{
                        //    Ch.Parent.Invoke(new VoidFunc(Ch.Series[0].Points.Clear), new object[] { });
                        //    foreach (PopulationElement p in Population.Population)
                        //    {
                        //        Ch.Parent.Invoke(new Func<double, double, int>(Ch.Series[0].Points.AddXY), new object[] { p.unfitness, p.fitness });
                        //    }
                        //}

                        if (ChartUnfitness != null)
                        {
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[0].Points.Clear), new object[] { });
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[1].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[0].Points.AddXY), new object[] { 0, new object[] { X[0].unfitness } });
                        }

                        if (ChartFitness != null)
                        {
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[0].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[1].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[0].Points.AddXY), new object[] { 0, new object[] { X[0].fitness } });
                        }
                        
                        for (int j = 0; j < dataColumnCount ; j++)
                        {
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Remove(Convert.ToString(j )); }));
                            data.Invoke(new Action(delegate { data.Columns.Remove(Convert.ToString(j )); }));
                        }
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                            data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                            data_vectors.Invoke(new Action(delegate { data_vectors.Rows.Add(ColCount - 1); }));
                            data.Invoke(new Action(delegate { data.Rows.Add(RowCount - 1); }));
                            for (int m = 0; m < ColCount; m++)
                                data_vectors[0, m].Value = X[0].Value[m];
                            sol = new int[RowCount];
                            for (int n = 0; n < RowCount; n++)
                            {
                                for (int m = 0; m < ColCount; m++)
                                    sol[n] += X[0].Value[m] * A[n][m];
                                data[0, n].Value = sol[n];
                            }
                        #endregion
                        int i = 1;
                        _min = Double.MaxValue;
                        do
                        {
                            EpsPopulationElement Y;

                            byte[] array = new byte[X.Count];
                            array = X[0].Value.ToArray();

                            bool found = false, feasible = false;
                            do
                            {
                                EpsPopulation = new EpsPopulationClass(this, eps, array);
                                Y = EpsPopulation.ElementAt(0);
                                for (int j = 1; j < EpsPopulation.Count; j++)
                                {
                                    do
                                    {
                                        if (MemoryList.ElementAt(MemoryList.Count-1).unfitness != EpsPopulation.ElementAt(j).unfitness || MemoryList.ElementAt(MemoryList.Count-1).fitness != EpsPopulation.ElementAt(j).fitness)
                                        {
                                            if (EpsPopulation.ElementAt(j).unfitness < Y.unfitness)
                                            {
                                                Y = EpsPopulation.ElementAt(j); feasible = true;
                                            }
                                            else if (EpsPopulation.ElementAt(j).unfitness == Y.unfitness)
                                            {
                                                if (EpsPopulation.ElementAt(j).fitness < Y.fitness)
                                                {
                                                    Y = EpsPopulation.ElementAt(j);
                                                    feasible = true;
                                                }
                                            }
                                            else feasible = true;
                                        }
                                        else { feasible = true; }
                                    } while (feasible != true);
                                }

                                    if (Y.unfitness < X[0].unfitness)
                                    {
                                        X[0] = X[0].Convert(Y, X[0]);
                                        MemoryList.Add(X[0]); eps = 4;
                                        found = true;
                                        break;
                                    }
                                    else if (Y.unfitness == X[0].unfitness)
                                    {
                                        if (Y.fitness < X[0].fitness)
                                        {
                                            X[0] = X[0].Convert(Y, X[0]);
                                            MemoryList.Add(X[0]); eps = 4;
                                            found = true;
                                            break;
                                        }
                                        else
                                        {
                                                eps *= 2;
                                        }
                                    }
                                    else
                                    {
                                            eps *= 2;
                                    }
                            } while (found != true);
                        

                            //  while (Population.Exists(X));

                            //#region Drawing
                            //if (Ch != null)
                            //{
                            //    Ch.Parent.Invoke(new VoidFunc(Ch.Series[1].Points.Clear), new object[] { });
                            //    Ch.Parent.Invoke(new Func<object, object[], int>(Ch.Series[1].Points.AddXY), new object[] { X.unfitness, new object[] { X.fitness } });
                            //}
                            //#endregion

                            if (SleepTime) Thread.Sleep(100);

                            //int index = Population.Replace(X);

                            #region Drawing
                            if (tbIteration != null)
                                tbIteration.Invoke(new VoidFunc<string>(text => tbIteration.Text = text), new object[] { String.Format("{0}/{1} Time: {2}", i, M, DateTime.Now.Subtract(T0)) });

                            if (ChartUnfitness != null)
                            {
                                //   List<PopulationElement> tmp = Population.Population.OrderBy(p => p.unfitness).ToList();
                                ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[0].Points.AddXY), new object[] { i, new object[] { X[0].unfitness } });
                                //    ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[1].Points.AddXY), new object[] { i, new object[] { tmp[tmp.Count - 1].unfitness } });
                            }

                            if (ChartFitness != null)
                            {
                                //    List<PopulationElement> tmp = Population.Population.OrderBy(p => p.fitness).ToList();
                                ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[0].Points.AddXY), new object[] { i, new object[] { X[0].fitness } });
                                //    ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[1].Points.AddXY), new object[] { i, new object[] { tmp[tmp.Count - 1].fitness } });
                            }
                            #endregion

                            if (X[0].unfitness == 0)
                            {
                                //_min = X[0].fitness;
                                // Console.Beep();
                                if (tbOptimum != null)
                                    tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = X[0].fitness.ToString()), new object[] { "" });

                            }
                            else
                                if (tbOptimum != null)
                                    tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = "No Solution found"), new object[] { "" });


                            //if (tbOptimum != null)
                            //    tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = _min.ToString()), new object[] { "" });

                            //}
                            //else
                            //if (_min == Double.MaxValue)
                            //{
                            //    if (tbOptimum != null)
                            //        tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = "No Solution found"), new object[] { "" });
                            //}
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                            data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                            for (int m = 0; m < ColCount; m++)
                                data_vectors[i , m].Value = X[0].Value[m];
                            sol = new int[RowCount];
                            for (int n = 0; n < RowCount; n++)
                            {
                                for (int m = 0; m < ColCount; m++)
                                    sol[n] += X[0].Value[m] * A[n][m];
                                data[i , n].Value = sol[n];
                            }
                            i++;
                        } while (i < M && Solving);
                        break;
                    }
                case AlgType.atHillClimbing:
                    {
                        PopulationClass X = new PopulationClass(this, 1);                        
                        #region Drawing
                        if (ChartUnfitness != null)
                        {
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[0].Points.Clear), new object[] { });
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[1].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[0].Points.AddXY), new object[] { 0, new object[] { X[0].unfitness } });
                        }

                        if (ChartFitness != null)
                        {
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[0].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[1].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[0].Points.AddXY), new object[] { 0, new object[] { X[0].fitness } });
                        }
                        {
                            for (int j = 0; j < dataColumnCount ; j++)
                            {
                                data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Remove(Convert.ToString(j )); }));
                                data.Invoke(new Action(delegate { data.Columns.Remove(Convert.ToString(j )); }));
                            }
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                            data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                            data_vectors.Invoke(new Action(delegate { data_vectors.Rows.Add(ColCount - 1); }));
                            data.Invoke(new Action(delegate { data.Rows.Add(RowCount - 1); }));
                            for (int m = 0; m < ColCount; m++)
                                data_vectors[0, m].Value = X[0].Value[m];
                            sol = new int[RowCount];
                            for (int n = 0; n < RowCount; n++)
                            {
                                for (int m = 0; m < ColCount; m++)
                                    sol[n] += X[0].Value[m] * A[n][m];
                                data[0, n].Value = sol[n];
                            }
                            
                        }
                        #endregion
                        int i = 1;
                        _min = Double.MaxValue;
                        do
                        {
                            EpsPopulationElement Y;
                            
                         //   X = Population.ElementAt((int)(new MyRandom()).Next(Population.Count));
                            //if (X.unfitness == 0 && _min == X.fitness)
                            //    X = Population.ElementAt((int)(new MyRandom()).Next(Population.Count));
                            //else
                            
                                byte[] array = new byte[X.Count];
                                array = X[0].Value.ToArray();

                                bool fEps = false;
                                do
                                {
                                    EpsPopulation = new EpsPopulationClass(this, eps, array);

                                    Y = EpsPopulation.ElementAt((int)((new MyRandom()).Next(EpsPopulation.Count)));
                                    bool found = false;
                                    do
                                    {
                                        if (Y.unfitness < X[0].unfitness)
                                        {
                                            X[0] = X[0].Convert(Y, X[0]); eps = 4;
                                            found = true;
                                            fEps = true; 
                                            break;
                                        }
                                        if (Y.unfitness == X[0].unfitness)
                                        {
                                            if (Y.fitness < X[0].fitness)
                                            {
                                                X[0] = X[0].Convert(Y, X[0]); eps = 4;
                                                found = true;
                                                fEps = true;
                                                break;
                                            }
                                            //else Y = EpsPopulation.ElementAt((int)(new MyRandom()).Next(EpsPopulation.Count));
                                            else
                                            {
                                                EpsPopulation.Delete(Y);
                                                if (EpsPopulation.Count != 0)
                                                {
                                                    Y = EpsPopulation.ElementAt((int)(new MyRandom()).Next(EpsPopulation.Count));
                                                }
                                                else 
                                                { 
                                                   eps *= 2; 
                                                   found = true; 
                                                }
                                            }
                                        }
                                        else
                                        {
                                            // Y = EpsPopulation.ElementAt((int)(new MyRandom()).Next(EpsPopulation.Count)); 
                                            EpsPopulation.Delete(Y);
                                            if (EpsPopulation.Count != 0)
                                            {
                                                Y = EpsPopulation.ElementAt((int)(new MyRandom()).Next(EpsPopulation.Count));
                                            }
                                            else
                                            { 
                                                eps *= 2;  
                                                found = true; 
                                            }
                                        }
                                    } while (found != true);
                                } while (fEps != true);
                            

                          //  while (Population.Exists(X));

                            //#region Drawing
                            //if (Ch != null)
                            //{
                            //    Ch.Parent.Invoke(new VoidFunc(Ch.Series[1].Points.Clear), new object[] { });
                            //    Ch.Parent.Invoke(new Func<object, object[], int>(Ch.Series[1].Points.AddXY), new object[] { X.unfitness, new object[] { X.fitness } });
                            //}
                            //#endregion

                            if (SleepTime) Thread.Sleep(100);

                            //int index = Population.Replace(X);

                            #region Drawing
                            if (tbIteration != null)
                                tbIteration.Invoke(new VoidFunc<string>(text => tbIteration.Text = text), new object[] { String.Format("{0}/{1} Time: {2}", i, M, DateTime.Now.Subtract(T0)) });

                            if (ChartUnfitness != null)
                            {
                            //   List<PopulationElement> tmp = Population.Population.OrderBy(p => p.unfitness).ToList();
                                ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[0].Points.AddXY), new object[] { i, new object[] {X[0].unfitness } });
                            //    ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[1].Points.AddXY), new object[] { i, new object[] { tmp[tmp.Count - 1].unfitness } });
                            }

                            if (ChartFitness != null)
                            {
                            //    List<PopulationElement> tmp = Population.Population.OrderBy(p => p.fitness).ToList();
                                ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[0].Points.AddXY), new object[] { i, new object[] { X[0].fitness } });
                            //    ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[1].Points.AddXY), new object[] { i, new object[] { tmp[tmp.Count - 1].fitness } });
                            }
                            #endregion
                            //PopulationElement[] PP = (from p in Population.Population where p.unfitness == 0 select p).ToArray();

                            //if (PP.Length > 0)
                            //{

                              //  PopulationElement Solvation = PP.OrderBy(p => p.fitness).ToList()[0];

                                sol = new int[RowCount];
                                for (int n = 0; n < RowCount; n++)
                                {
                                    for (int m = 0; m < ColCount; m++)

                                        sol[n] += X[0].Value[m] * A[n][m];
                                    data[data.ColumnCount - 1, n].Value = sol[n];
                                }

                                if (X[0].unfitness == 0)
                                {
                                    //_min = X[0].fitness;
                                   // Console.Beep();
                                    if (tbOptimum != null)
                                        tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = X[0].fitness.ToString()), new object[] { "" });

                                }
                                else
                                if (tbOptimum != null)
                                    tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = "No Solution found"), new object[] { "" });
                                

                                //if (tbOptimum != null)
                                //    tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = _min.ToString()), new object[] { "" });

                            //}
                            //else
                                //if (_min == Double.MaxValue)
                                //{
                                //    if (tbOptimum != null)
                                //        tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = "No Solution found"), new object[] { "" });
                                //}
                                //data_vectors.Columns.Add(Convert.ToString(i + 1), Convert.ToString(i + 1));
                                data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                                data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                                for (int m = 0; m < ColCount; m++)
                                    data_vectors[i , m].Value = X[0].Value[m];
                                sol = new int[RowCount];
                                for (int n = 0; n < RowCount; n++)
                                {
                                    for (int m = 0; m < ColCount; m++)
                                        sol[n] += X[0].Value[m] * A[n][m];
                                    data[i, n].Value = sol[n];
                                }
                            i++;
                        } while (i < M && Solving);
                            break;}
                    
                    case AlgType.atGeneticAlgorithm:
                    {
                        Population = new PopulationClass(this, 100);

                        #region Drawing
                        if (Ch != null)
                        {
                            Ch.Parent.Invoke(new VoidFunc(Ch.Series[0].Points.Clear), new object[] { });
                            foreach (PopulationElement p in Population.Population)
                            {
                                Ch.Parent.Invoke(new Func<double, double, int>(Ch.Series[0].Points.AddXY), new object[] { p.unfitness, p.fitness });
                            }
                        }

                        if (ChartUnfitness != null)
                        {
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[0].Points.Clear), new object[] { });
                            ChartUnfitness.Invoke(new VoidFunc(ChartUnfitness.Series[1].Points.Clear), new object[] { });
                        }

                        if (ChartFitness != null)
                        {
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[0].Points.Clear), new object[] { });
                            ChartFitness.Invoke(new VoidFunc(ChartFitness.Series[1].Points.Clear), new object[] { });
                        }
                        for (int j = 0; j < dataColumnCount; j++)
                        {
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Remove(Convert.ToString(j)); }));
                            data.Invoke(new Action(delegate { data.Columns.Remove(Convert.ToString(j)); }));
                        }
                        data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                        data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(0), "Iteration № 0"); }));
                        data_vectors.Invoke(new Action(delegate { data_vectors.Rows.Add(ColCount - 1); }));
                        data.Invoke(new Action(delegate { data.Rows.Add(RowCount - 1); }));
                        #endregion
                        int i = 0;
                        _min = Double.MaxValue;
                        do
                        {
                            PopulationElement[] Parents = Population.GetParents();
                            PopulationElement Child;

                            do
                            {
                                /*if ((from p in Population where p.unfitness == 0 select p.unfitness).Count() == Population.Count)
                                {
                                    Child = Parents[(int)(new MyRandom()).Next(2)];
                                }
                                else*/
                                Child = PopulationClass.Crossover(Parents[0], Parents[1]);
                                PopulationClass.Mutation(ref Child, DynamicMutation);
                                if (UsingHeuristicOperator) PopulationClass.HeuristicFeasibilityOperator(ref Child);
                            } while (Population.Exists(Child) /*|| Child.unfitness > Population.Population.OrderByDescending(p => p.unfitness).ToList()[0].unfitness*/);

                            #region Drawing
                            if (Ch != null)
                            {
                                Ch.Parent.Invoke(new VoidFunc(Ch.Series[1].Points.Clear), new object[] { });
                                Ch.Parent.Invoke(new Func<object, object[], int>(Ch.Series[1].Points.AddXY), new object[] { Child.unfitness, new object[] { Child.fitness } });
                            }
                            #endregion

                            if (SleepTime) Thread.Sleep(100);

                            int index = Population.Replace(Child);

                            #region Drawing
                            if (Ch != null && index > -1)
                            {
                                Ch.Parent.Invoke(new VoidFunc<int>(Ch.Series[0].Points.RemoveAt), new object[] { index });
                                Ch.Parent.Invoke(new VoidFunc<int, object, object[]>(Ch.Series[0].Points.InsertXY), new object[] { index, Child.unfitness, new object[] { Child.fitness } });
                            }

                            if (tbIteration != null)
                                tbIteration.Invoke(new VoidFunc<string>(text => tbIteration.Text = text), new object[] { String.Format("{0}/{1} Time: {2}", i, M, DateTime.Now.Subtract(T0)) });

                            if (ChartUnfitness != null)
                            {
                                List<PopulationElement> tmp = Population.Population.OrderBy(p => p.unfitness).ToList();
                                ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[0].Points.AddXY), new object[] { i, new object[] { tmp[0].unfitness } });
                                ChartUnfitness.Invoke(new Func<object, object[], int>(ChartUnfitness.Series[1].Points.AddXY), new object[] { i, new object[] { tmp[tmp.Count - 1].unfitness } });
                            }

                            if (ChartFitness != null)
                            {
                                List<PopulationElement> tmp = Population.Population.OrderBy(p => p.fitness).ToList();
                                ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[0].Points.AddXY), new object[] { i, new object[] { tmp[0].fitness } });
                                ChartFitness.Invoke(new Func<object, object[], int>(ChartFitness.Series[1].Points.AddXY), new object[] { i, new object[] { tmp[tmp.Count - 1].fitness } });
                            }
                            sol = new int[RowCount];
                            //for (int n = 0; n < RowCount; n++)
                            //{
                            //    for (int m = 0; m < ColCount; m++)

                            //        sol[n] += Solvation.Value[m] * A[n][m];
                            //    data[data.ColumnCount - 1, n].Value = sol[n];
                            //}
                            List<PopulationElement> _tmp = Population.Population.OrderBy(p => p.unfitness).ToList();
                            data_vectors.Invoke(new Action(delegate { data_vectors.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                            data.Invoke(new Action(delegate { data.Columns.Add(Convert.ToString(i), "Iteration №" + i); }));
                            for (int m = 0; m < ColCount; m++)
                                data_vectors[i, m].Value = _tmp[0].Value[m];
                            sol = new int[RowCount];
                            for (int n = 0; n < RowCount; n++)
                            {
                                for (int m = 0; m < ColCount; m++)
                                    sol[n] += _tmp[0].Value[m] * A[n][m];
                                data[i, n].Value = sol[n];
                            }
                            #endregion
                            PopulationElement[] PP = (from p in Population.Population where p.unfitness == 0 select p).ToArray();
                                                       
                            if (PP.Length > 0)
                            {

                                PopulationElement Solvation = PP.OrderBy(p => p.fitness).ToList()[0];

                                //for (int n = 0; n < Solvation.Count; n++)
                                //    data[data.ColumnCount - 1, n].Value = Solvation.Value[n];

                                
                                if (Solvation.fitness < _min)
                                {
                                    _min = Solvation.fitness;
                                   // Console.Beep();
                                }

                                if (tbOptimum != null)
                                    tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = _min.ToString()), new object[] { "" });

                            }
                            else
                                if (_min == Double.MaxValue)
                                {
                                    if (tbOptimum != null)
                                        tbOptimum.Invoke(new VoidFunc<string>(text => tbOptimum.Text = "No Solution found"), new object[] { "" });
                                }
                            i++;
                        } while (i < M && Solving/* && (from p in Population where p.unfitness==0 select p.unfitness).Count()!=Population.Count*/);
                        break;
                    }
            }
        }

        public static double Fitness(byte[] x, byte[][] A, int[] c, out double fitness, out double unfitness)
        {
            fitness = 0;
            unfitness = 0;
            int _ColCount = x.Length;
            int _RowCount = A.Length;
            int i, j;
            for (j = 0; j < _ColCount; j++)
            {
                fitness += c[j] * x[j];
            }

            for (i = 0; i < _RowCount; i++)
            {
                //int n = A[i].Length;
                double tmp = 0;
                for (j = 0; j < _ColCount; j++)
                {
                    tmp += x[j] * A[i][j];
                }

                unfitness += Math.Abs(tmp - 1);
            }

            return fitness + unfitness;
        }
        #endregion
    }
}
