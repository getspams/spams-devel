
/* Software SPAMS v2.1 - Copyright 2009-2011 Julien Mairal 
 *
 * This file is part of SPAMS.
 *
 * SPAMS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SPAMS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SPAMS.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PROJECT_H
#define PROJECT_H

#include <linalg.h>
#include <limits>

#define EPSILON_MAXFLOW 1e-10

//#define VERBB
//#define VERB2


int num_relabels;
int num_pushes;
int num_global_relabels;
int num_gap_relabels;
bool global_heuristic = true;
bool gap_heuristic = true;
bool cap_heuristic = true;
bool price_heuristic = true;
bool price_refine_heuristic = false;

//typedef std::list<int> list_int;
//typedef std::list<int>::const_iterator const_iterator_int;
#include <list.h>
typedef List<int> list_int;
typedef ListIterator<int> const_iterator_int;

Timer tglobal1, tglobal2, tglobal3;

template <typename T>
bool compare_abs (T first, T second) {
   return abs<T>(first) >= abs<T>(second);
}
template <typename T>
T inline project_tree_l1(T* variables, const int n, const T lambda);

template <typename T>
T inline project_tree_l1(T* X, const int n, const T lambda) {
   if (lambda==0) return INFINITY;
   T* prU = X;
   T sum = 0;
   int sum_card = n;
   for (int i = 0; i<sum_card; ++i) {
      if (X[i]) {
         sum += X[i];
      } else {
         swap(X[i],X[--sum_card]);
         --i;
      }
   }
   if (sum < lambda) {
      memset(X,0,sum_card*sizeof(T));
      return 0;
   }
   int sizeU = sum_card;
   sum_card = 0;
   sum=0;

   while (sizeU > 0) {
      // put the pivot in prU[0]
      swap(prU[0],prU[sizeU/2]);
      int sizeG=1;
      T sumG=prU[0];

      for (int i = 1; i<sizeU; ++i) {
         if (prU[i] >= prU[0]) {
            sumG += prU[i];
            swap(prU[sizeG++],prU[i]);
         }
      }

      T new_sum=sum+sumG;
      int new_card=sum_card+sizeG;
      if (new_sum - prU[0]*new_card <= lambda) {
         sum_card = new_card;
         sum = new_sum;
         prU +=sizeG;
         sizeU -= sizeG;
      } else {
         ++prU;
         sizeU = sizeG-1;
      }
   }
   T thrs = MAX(0,(sum-lambda)/sum_card);
   for (int i = 0; i<n; ++i) 
      X[i] = MIN(X[i],thrs);
   return thrs;
};

template <typename T> class Tree_Seq {

   public:
      Tree_Seq();
      ~Tree_Seq();

      void inline create_tree(const int N_variables, int* own_variables,
            int* N_own_variables, T* lambda, mwSize* groups_ir, mwSize* groups_jc,
            const int N_groups, const int root_node = 0);

      int inline perform_order(const int current_node, const int pointer);
      int inline perform_dfs(const int current_node, const int pointer);

      void inline proj(Vector<T>& input, const bool l1 = false,
            const T fact = 1.0);
      void inline proj_zero(Vector<T>& input, const T fact = 1.0);

      void inline proj_weighted_linf(Vector<T>& input, const Vector<T>& weights, const T fact = 1.0);

      T inline val_norm(const T* pr_alpha, const int current_node, const bool l1 = false);
      T inline val_norm2(const T* pr_alpha, const int current_node, T& tmp, const bool l1 = false);
      T inline val_zero(const T* pr_alpha, const int current_node);
      T inline val_zero2(const T* pr_alpha, const int current_node, bool& tmp);
      T inline dual_norm_inf(const Vector<T>& input);
      void inline sub_grad(const Vector<T>& input,  Vector<T>& output, const bool linf);

   private:
      int _N_groups;
      int _N_variables;
      T* _lambda;
      T* _thrs;
      T* _variables;
      T* _work;
      int* _size_variables;
      int* _pr_variables;
      int* _size_own_variables;
      int* _pr_own_variables;
      int* _order;
      int* _order_dfs;
      mwSize* _groups_ir;
      mwSize* _groups_jc;

};

template <typename T> 
Tree_Seq<T>::Tree_Seq() {
   _lambda=NULL;
   _thrs=NULL;
   _work=NULL;
   _variables=NULL;
   _N_groups=0;
   _N_variables=0;
   _size_variables=NULL;
   _pr_variables=NULL;
   _size_own_variables=NULL;
   _order=NULL;
   _order_dfs=NULL;
   _groups_ir=NULL;
   _groups_jc=NULL;
};

template <typename T>
Tree_Seq<T>::~Tree_Seq() {
   delete[](_thrs);
   delete[](_work);
   delete[](_variables);
   delete[](_size_variables);
   delete[](_pr_variables);
   delete[](_order);
   delete[](_order_dfs);
};

template <typename T>
void inline Tree_Seq<T>::create_tree(const int N_variables, int* own_variables,
      int* N_own_variables, T* lambda, mwSize* groups_ir,mwSize* groups_jc,
      const int N_groups, const int root_node) {
   _N_groups=N_groups;
   _N_variables=N_variables;
   _lambda=lambda;
   _thrs=new T[N_groups];
   _variables=new T[N_variables];
   _size_variables=new int[N_groups];
   _pr_variables=new int[N_groups];
   _size_own_variables=N_own_variables;
   _pr_own_variables=own_variables;
   _order=new int[N_groups];
   _order_dfs=new int[N_groups];
   _groups_ir=groups_ir;
   _groups_jc=groups_jc;
   this->perform_order(root_node,0);
   this->perform_dfs(root_node,0);
   _work = new T[MAX(N_groups,N_variables)];
};

template <typename T>
int inline Tree_Seq<T>::perform_order(const int current_node, const int pointer) {
   int cur_pointer=pointer;
   _size_variables[current_node]=_size_own_variables[current_node];
   _pr_variables[current_node]=_pr_own_variables[current_node];
   for (int i = _groups_jc[current_node];  i<_groups_jc[current_node+1]; ++i) {
      cur_pointer=this->perform_order(_groups_ir[i],cur_pointer);
      _size_variables[current_node]+=_size_variables[_groups_ir[i]];
      _pr_variables[current_node]= MIN(_pr_variables[current_node],_pr_variables[_groups_ir[i]]);
   }
   _order[cur_pointer]=current_node;
   return cur_pointer+1;
};

template <typename T>
int inline Tree_Seq<T>::perform_dfs(const int current_node, const int pointer) {
   int cur_pointer=pointer;
   _order_dfs[cur_pointer++]=current_node;
   for (int i = _groups_jc[current_node];  i<_groups_jc[current_node+1]; ++i) {
      cur_pointer=this->perform_dfs(_groups_ir[i],cur_pointer);
   }
   return cur_pointer;
};

// could be faster
template <typename T>
T inline Tree_Seq<T>::val_norm(const T* pr_alpha, const int current_node, const bool l1) {
   T tmp=0;
   return this->val_norm2(pr_alpha,current_node,tmp,l1);
};

// fast version
template <typename T>
T inline Tree_Seq<T>::val_norm2(const T* pr_alpha, const int current_node, T& tmp, const bool l1) {
   T sum=0;
   for (int i = _groups_jc[current_node];  i<_groups_jc[current_node+1]; ++i) {
      T tmp2=0;
      sum+=this->val_norm2(pr_alpha,_groups_ir[i],tmp2,l1);
      tmp= l1 ? MAX(tmp,tmp2) : tmp+tmp2;
   }
   if (l1) {
      for (int i = 0; i<_size_own_variables[current_node]; ++i)
         tmp=MAX(abs<T>(pr_alpha[i+_pr_variables[current_node]]),tmp);
      sum+=_lambda[current_node]*tmp;
   } else {
      tmp += cblas_dot<T>(_size_own_variables[current_node],const_cast<T*>(pr_alpha+_pr_variables[current_node]),1,const_cast<T*>(pr_alpha+_pr_variables[current_node]),1);
      //tmp += cblas_dot<T>(_size_own_variables[current_node],pr_alpha+_pr_variables[current_node],1,pr_alpha+_pr_variables[current_node],1);
      sum+=_lambda[current_node]*sqrt(tmp);
   }
   return sum;
};

template <typename T>
T inline Tree_Seq<T>::val_zero2(const T* pr_alpha, const int current_node, bool& tmp) {
   T sum=0;
   for (int i = _groups_jc[current_node];  i<_groups_jc[current_node+1]; ++i) {
      bool tmp2=false;
      sum+=this->val_zero2(pr_alpha,_groups_ir[i],tmp2);
      tmp = tmp || tmp2;
   }
   for (int i = 0; i<_size_own_variables[current_node]; ++i)
      tmp= (tmp || pr_alpha[i+_pr_variables[current_node]] != 0);
   if (tmp)
      sum+=_lambda[current_node];
   return sum;
};


template <typename T>
T inline Tree_Seq<T>::val_zero(const T* pr_alpha, const int current_node) {
   bool tmp = false;
   return this->val_zero2(pr_alpha,current_node,tmp);
};

template <typename T>
void inline Tree_Seq<T>::sub_grad(const Vector<T>& input, Vector<T>& output, const bool linf) {
   output.setZeros();
   if (linf) {
      for (int i = 0; i<_N_groups; ++i) {
         const T* pr = input.rawX()+_pr_variables[i];
         int imax=cblas_iamax<T>(_size_variables[i],const_cast<T*>(pr),1);
      //   int imax=cblas_iamax<T>(_size_variables[i],pr,1);
         T max=pr[imax];
         int num_max=0;
         for (int j = 0; j<_size_variables[i];++j) {
            if (abs<T>(max-abs<T>(pr[j])) < 1e-10) ++num_max;
         }
         T add=T(1.0)/num_max;
         for (int j = 0; j<_size_variables[i];++j) {
            if (abs<T>(max-abs<T>(pr[j])) < 1e-10 && input[_pr_variables[i]+j]) {
               output[_pr_variables[i]+j]+= input[_pr_variables[i]+j] > 0 ? add : -add;
            }
         }
      }
   } else {
      for (int i = 0; i<_N_groups; ++i) {
         T nrm=cblas_nrm2<T>(_size_variables[i],input.rawX()+_pr_variables[i],1);
         if (nrm > 0) {
            cblas_axpy<T>(_size_variables[i],T(1.0)/nrm,input.rawX()+_pr_variables[i],1,output.rawX()+_pr_variables[i],1);
//         } else {
//            T num=T(1.0)/sqrt(static_cast<T>(_size_variables[i]));
//            for (int j = 0; j<_size_variables[i]; ++j) {
//               output[_pr_variables[i]+j]+=num;
//            }
         }
      }
   }
};

template <typename T>
void inline Tree_Seq<T>::proj_weighted_linf(Vector<T>& input, const Vector<T>& weights, 
      const T fact) {
   Vector<T> weig;
   weig.copy(weights);
   weig.inv();
   cblas_copy<T>(input.n(),input.rawX(),1,_variables,1);
   Vector<T> tmp, tmpw;
   for (int i = 0; i<_N_groups; ++i) {
      const int node=_order[i];
      Vector<T> out;
      tmp.setData(_variables+_pr_variables[node],_size_variables[node]);
      tmpw.setData(weig.rawX()+_pr_variables[node],_size_variables[node]);
      tmp.l1project_weighted(out,tmpw,fact,true);
      cblas_copy<T>(out.n(),out.rawX(),1,_variables+_pr_variables[node],1);
      tmp.copy(out);
   }
   cblas_copy<T>(input.n(),_variables,1,input.rawX(),1);
};

template <typename T>
void inline Tree_Seq<T>::proj(Vector<T>& input, const bool l1,
      const T fact) {
   if (l1) {
      vAbs<T>(input.n(),input.rawX(),_variables);
      for (int i = 0; i<_N_groups; ++i) {
         const int node=_order[i];
         _thrs[node] = project_tree_l1(_variables+_pr_variables[node],_size_variables[node],
               _lambda[node]*fact);
      }
      cblas_copy<T>(input.n(),input.rawX(),1,_variables,1);
      for (int i = 0; i<_N_groups; ++i) {
         const int node=_order_dfs[i];
         if (_thrs[node] == 0) {
            memset(_variables+_pr_own_variables[node],0,_size_own_variables[node]*sizeof(T));
            for (int j = _groups_jc[node];  j<_groups_jc[node+1]; ++j) {
               _thrs[_groups_ir[j]]=0;
            }
         } else {
            for (int j = 0; j<_size_own_variables[node]; ++j) {
               T tmp = _variables[_pr_own_variables[node]+j];
               _variables[_pr_own_variables[node]+j] = tmp > _thrs[node] ? _thrs[node] :
                  tmp < -_thrs[node] ? -_thrs[node] : tmp;
            }
            for (int j = _groups_jc[node];  j<_groups_jc[node+1]; ++j) {
               _thrs[_groups_ir[j]]= MIN(_thrs[_groups_ir[j]],_thrs[node]);
            }
         }
      }
   } else {
      cblas_copy<T>(input.n(),input.rawX(),1,_variables,1);
      for (int i = 0; i<_N_groups; ++i) {
         const int node=_order[i];
         _work[node]=0;
         for (int j = 0; j<_size_own_variables[node]; ++j)
            _work[node]+=_variables[_pr_own_variables[node]+j]*_variables[_pr_own_variables[node]+j];
         for (int j = _groups_jc[node]; j<_groups_jc[node+1];++j)
            _work[node] += _work[_groups_ir[j]];
         _thrs[node] = MAX(0,1-fact*_lambda[node]/sqrt(_work[node]));
         _work[node]*= _thrs[node]*_thrs[node];
      }
      for (int i = 0; i<_N_groups; ++i) {
         const int node=_order_dfs[i];
         if (_thrs[node] == 0) {
            memset(_variables+_pr_own_variables[node],0,_size_own_variables[node]*sizeof(T));
            for (int j = _groups_jc[node];  j<_groups_jc[node+1]; ++j) {
               _thrs[_groups_ir[j]]=0;
            }
         } else {
            for (int j = 0; j<_size_own_variables[node]; ++j) 
               _variables[_pr_own_variables[node]+j] *= _thrs[node];
            for (int j = _groups_jc[node];  j<_groups_jc[node+1]; ++j) {
               _thrs[_groups_ir[j]] *= _thrs[node];
            }
         }
      }
   }
   cblas_copy<T>(input.n(),_variables,1,input.rawX(),1);
};

template <typename T>
void inline Tree_Seq<T>::proj_zero(Vector<T>& input, const T fact) {
   cblas_copy<T>(input.n(),input.rawX(),1,_variables,1);
   for (int i = 0; i<_N_groups; ++i) {
      const int node=_order[i];
      _work[node]=0;
      for (int j = 0; j<_size_own_variables[node]; ++j)
         _work[node]+=_variables[_pr_own_variables[node]+j]*_variables[_pr_own_variables[node]+j];
      _work[node] *= -0.5;
      _work[node] += fact*_lambda[node];
      for (int j = _groups_jc[node]; j<_groups_jc[node+1];++j)
         _work[node] += _work[_groups_ir[j]];
      if (_work[node] > 0) _work[node]=0;
   }
   for (int i = 0; i<_N_groups; ++i) {
      const int node=_order_dfs[i];
      if (_work[node] == 0) {
         memset(_variables+_pr_own_variables[node],0,_size_own_variables[node]*sizeof(T));
         for (int j = _groups_jc[node];  j<_groups_jc[node+1]; ++j) {
            _work[_groups_ir[j]]=0;
         }
      } 
   }
   cblas_copy<T>(input.n(),_variables,1,input.rawX(),1);
}

template <typename T>
T inline Tree_Seq<T>::dual_norm_inf(const Vector<T>& input) {
   tglobal1.reset();
   tglobal2.reset();
   tglobal3.reset();
   for (int i = 0; i<_N_groups; ++i) {
      _thrs[_order[i]]=INFINITY;
   }
   T tau=0;
   T sum_variables=INFINITY;
   T total=input.asum();
   while (_thrs[0] > EPSILON) {
      T old_thrs=_thrs[0];
      vAbs<T>(_N_variables,input.rawX(),_variables);
      list_int nodes;
      nodes.push_front(0);
      list_int ordered_nodes;
      T sum_weights=0;
      sum_variables=total;
      while (!nodes.empty()) {
         const int node=nodes.front();
         nodes.pop_front();
         sum_weights+=_lambda[node];
         for (int j = _groups_jc[node];  j<_groups_jc[node+1]; ++j) 
            if (_thrs[_groups_ir[j]] > EPSILON) {
               nodes.push_front(_groups_ir[j]);
            } else {
               sum_variables -= cblas_asum<T>(_size_variables[_groups_ir[j]],_variables+_pr_variables[_groups_ir[j]],1);
               memset(_variables+_pr_variables[_groups_ir[j]],0,_size_variables[_groups_ir[j]]*sizeof(T));
            }
         ordered_nodes.push_front(node);
      }
      tau=sum_variables/sum_weights;
      for (const_iterator_int it = ordered_nodes.begin(); it != ordered_nodes.end(); ++it) {
         const int node=*it;
         _thrs[node] = project_tree_l1(_variables+_pr_variables[node],_size_variables[node],_lambda[node]*tau);
      }
      if (_thrs[0] >= old_thrs) break;
   }
   return tau;
};


template <typename T> class MaxFlow {

   public:
      MaxFlow(const int N, const int* num_edges, const int s, const int t);
      ~MaxFlow();

      void inline add_edge(const int u, const int v, const T cu, const T cv);

      void inline discharge(const list_int& component, const int u,const int max_label);

      void inline global_relabelling();

      void inline perform_maxflow();

      void inline print_excess();
      void inline print_labels();

      void inline deactivate();
      void inline deactivate(const list_int& component);

      T inline getMaxFlow() const { return _excess[_t]; };

      void inline extractConnexComponents(std::list< list_int* >& connex_components);
      T inline project(const list_int& component, 
            const T* variables_in,T* variables_out,T* work,
            const int Ng);
      T inline project_weighted(const list_int& component, 
            const T* variables_in,T* variables_out,T* work, const T* weights,
            const int Ng);
      T inline project_box(const list_int& component, const T* variables_in,T*
            variables_out,T* work, bool& fusion, const int Ng);

      T inline flow_component(const list_int& component, const int Ng);
      bool inline splitComponent(const list_int& component,
            std::list< list_int* >& connex_components, const int Ng, bool* positive,const bool addpos = true);
      void inline reset_component(const list_int& component);
      void inline perform_maxflow_component(const list_int& component);
      void inline component_relabelling(const list_int& component,
            const int max_label, const bool force);
      void inline gap_relabelling(const list_int& component, const int gap, const int max_label);
      void inline component_gap(const list_int& component);
      void inline update_capacities(const list_int& component,T* work);
      void inline set_capacities_variables(const T* cap,const int Nv, const int Ng);
      void inline set_capacities_groups(const list_int& component,
            const Vector<T>& weights,const T lambda, const int Ng);
      void inline update_capacities_aux(const int node,T* work);
      void inline restore_capacities(const list_int& component);
      T inline compute_thrs_project_l1(T* X, const int n, const T lambda);

      bool inline check_flow();
      void inline restore_capacities();
      void inline restore_flow();
      void inline reset_flow();
      void inline scale_flow(const T scal);
      void inline save_capacities();
      void inline save_flow();
      void inline set_weights(const T lambda);
      void inline set_weights(const T* weights, const T lambda);
      void inline print_graph();
      inline void init_split_variables(SpMatrix<T>& splitted_w, const int Ng, const int Nv);
      inline void init_split_variables_aux(const int node, int& current_counter, Vector<int>& counter_node, list_int** splitted_w,
            const int Ng, const int Nv);
      void inline print_component(const list_int& component);
      void inline print_component2(const list_int& component);
      void inline print_sink();
      void inline print_graph_aux(const int node);
      T inline norm(const T* variables, T* work, const T* weights, const int Ng, const bool linf = true);
      inline int nzmax() const { return _nzmax; };
      void inline sub_gradient(const Vector<T>& input, Vector<T>& output, const Vector<T>& weights, const int Ng); 
      void inline sub_gradient_aux(const Vector<T>& input, Vector<T>& output, const Vector<T>& weights,
            const int node, list_int& list, const int Ng); 

   private:
      int _N;
      int _s;
      int _t;

      int* _labels;
      T* _excess;
      T* _copyexcess;
      bool* _seen;
      bool* _active;

      int* _max_num_edges;
      int* _current_edges;
      int* _num_edges;
      int* _pr_node;
      int _nzmax;

      int* _children;
      int* _reverse_address;
      T* _capacity;
      T* _copycapacity;
      T* _flow;
      T* _copyflow;

      int _current_max_label;
      list_int** _active_nodes;
      int* _all_nodes;
};

template <typename T>
MaxFlow<T>::MaxFlow(const int N, const int* num_edges, const int s, const int t) {
   _N=N;
   _s=s;
   _t=t;

   _labels=new int[N];
   memset(_labels,0,N*sizeof(int));
   _excess=new T[N];
   memset(_excess,0,N*sizeof(T));
   _excess[_s]=INFINITY;
   _seen=new bool[N];
   _active=new bool[N];
   _num_edges=new int[N];
   _current_edges=new int[N];
   memset(_num_edges,0,N*sizeof(int));
   memset(_current_edges,0,N*sizeof(int));
   _max_num_edges=new int[N];
   for (int i = 0; i<N; ++i) _max_num_edges[i]=num_edges[i];
   _pr_node=new int[N+1];
   _pr_node[0]=0;
   for (int i = 1; i<=N; ++i) _pr_node[i]=_pr_node[i-1]+_max_num_edges[i-1];
   _nzmax=_pr_node[N];
   _children= new int[_nzmax];
   _reverse_address= new int[_nzmax];
   _capacity= new T[_nzmax];
   _copycapacity= new T[_nzmax];
   _flow= new T[_nzmax];
   memset(_flow,0,_nzmax*sizeof(T));
   _current_max_label=0;
   _active_nodes = new list_int*[N+1];
   _all_nodes= new int[N+1];
   for (int i = 0; i<=N; ++i) _active_nodes[i]=new list_int();
};

template <typename T>
MaxFlow<T>::~MaxFlow() {
   delete[](_labels);
   delete[](_excess);
   delete[](_seen);
   delete[](_active);
   delete[](_num_edges);
   delete[](_current_edges);
   delete[](_max_num_edges);
   delete[](_children);
   delete[](_reverse_address);
   delete[](_capacity);
   delete[](_copycapacity);
   delete[](_flow);
   for (int i = 0; i<=_N; ++i) delete(_active_nodes[i]);
   delete[](_active_nodes);
   delete[](_all_nodes);
   delete[](_pr_node);
};

template <typename T>
void inline MaxFlow<T>::add_edge(const int u, const int v, const T Cu, const T Cv) {
   if (u != v) {
      const int pu=_pr_node[u];
      const int pv=_pr_node[v];
      const int nu=_num_edges[u]+pu;
      const int nv=_num_edges[v]+pv;
      _children[nu]=v;
      _children[nv]=u;
      _capacity[nu]=Cu;
      _capacity[nv]=Cv;
      _reverse_address[nu]=nv;
      _reverse_address[nv]=nu;
      _num_edges[u]++;
      _num_edges[v]++;
   }
};

template <typename T>
void inline MaxFlow<T>::discharge(const list_int& component, const int u, const int max_label) {
#ifdef VERBB
   cerr << "Discharge " << u << endl;
#endif
   const T* capacity=_capacity+_pr_node[u];
   T* flow=_flow+_pr_node[u];
   const int* children=_children+_pr_node[u];
   const int* reverse=_reverse_address+_pr_node[u];
   int pr=0;
   const int curr=_current_edges[u];
   const int nn=_num_edges[u];

   int m = max_label;
   while (_excess[u] > EPSILON_MAXFLOW && pr < nn) {
      const int num_c=(pr+curr) % nn;
      const int v=children[num_c];
      if (capacity[num_c] > flow[num_c]) {
         if (_labels[u] > _labels[v]) {
            // push
            const T delta=MIN(_excess[u],capacity[num_c]-flow[num_c]);
            _excess[u]-=delta;
            flow[num_c]+=delta;
#ifdef VERBB
            cerr << "Send " << delta << " from " << u << " to " << v << endl;
#endif
            /// add v to the list of active nodes
            if (!_active[v] && v != _t) {
               _active_nodes[_labels[v]]->push_back(v);
               _active[v]=true;
            }
            _excess[v]+=delta;
            _flow[reverse[num_c]]-=delta;
         } else {
            m=MIN(m,_labels[v]+1);
         }
      }
      ++pr;
   }
   num_relabels++;
   if (_excess[u]  > EPSILON_MAXFLOW) {
      /// relabel: 
      if (gap_heuristic) {
         _all_nodes[_labels[u]]--;
         if (_all_nodes[_labels[u]]==0) {
            this->gap_relabelling(component,_labels[u],max_label);
            _labels[u]=max_label;
         } else { 
#ifdef VERBB
            cerr << "relabel " << u << " from " << _labels[u] << " to " << MIN(m,max_label) << endl;
#endif
            _labels[u]=MIN(m,max_label);
            _all_nodes[_labels[u]]++;         
         }
      } else {
#ifdef VERBB
         cerr << "relabel " << u << " from " << _labels[u] << " to " << MIN(m,max_label) << endl;
#endif
         _labels[u]=MIN(m,max_label);
      }
   } else {
      _excess[u]=0;
      _current_edges[u]=((pr+curr) % nn);
   }
};

template <typename T>
void inline MaxFlow<T>::perform_maxflow() {
   int counter=1;
   while (_current_max_label > 0 || !_active_nodes[0]->empty()) {
      if (counter % 2*_N == 0) this->global_relabelling();
      if (_active_nodes[_current_max_label]->empty()) {
         _current_max_label--;
      } else {
         const int current_node=_active_nodes[_current_max_label]->front();
         _active_nodes[_current_max_label]->pop_front();
         _active[current_node]=false;
         if (_excess[current_node] > EPSILON_MAXFLOW) {
            this->discharge(NULL,current_node,_N);
            if (_excess[current_node] > EPSILON_MAXFLOW && _labels[current_node] < _N) {
               _active_nodes[_labels[current_node]]->push_back(current_node);
               _active[current_node]=true;
               if (_labels[current_node] > _current_max_label) 
                  _current_max_label=_labels[current_node];
            }
         } else {
            _excess[current_node]=0;
         }
      }
      ++counter;
   }
}

template <typename T>
void inline MaxFlow<T>::print_excess() {
   cerr << "Excess: " << endl;
   for (int i= 0; i<_N; ++i) {
      cerr << _excess[i] << " ";
   }
   cerr << endl;

};


template <typename T>
void inline MaxFlow<T>::deactivate() {
   for (int i= 0; i<_N; ++i) {
      _seen[i]=true;
      _active[i]=false;
      _labels[i]=_N;
   }
};

template <typename T>
void inline MaxFlow<T>::deactivate(const list_int& component) {
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      _seen[*it]=true;
      _active[*it]=false;
      _labels[*it]=_N;
   }
};


template <typename T>
void inline MaxFlow<T>::print_labels() {
   cerr << "Labels: " << endl;
   for (int i= 0; i<_N; ++i)
      cerr << _labels[i] << " ";
   cerr << endl;

};

template <typename T>
void inline MaxFlow<T>::print_graph() {
   cerr << "Number of nodes: " << _N << endl;
   cerr << "Source: " << _s << endl;
   cerr << "Sink: " << _t << endl;
   for (int i = 0; i<_N; ++i) _seen[i]=false;
   this->print_graph_aux(_s);
};

template <typename T>
void inline MaxFlow<T>::init_split_variables(SpMatrix<T>& splitted_w, const int Ng, const int Nv) {
   for (int i = 0; i<_N; ++i) _seen[i]=false;
   Vector<int> count(Ng);
   int current = 0;
   list_int** tab_list = new list_int*[Ng];
   for (int i = 0; i<Ng; ++i) tab_list[i] = new list_int();
   this->init_split_variables_aux(_s,current,count,tab_list,Ng,Nv);
   int nzmax = 0;
   for (int i = 0; i<Ng; ++i) {
      nzmax += tab_list[i]->size();      
   }
   /// assumes memory is not an issue
   splitted_w.resize(Nv,Ng,nzmax);
   int* pB = splitted_w.pB();
   int* r = splitted_w.r();
   T* v = splitted_w.v();
   pB[0]=0;
   int counter=0;
   for (int i = 0; i<Ng; ++i) {
      pB[i+1]=pB[i]+tab_list[i]->size();
      for (const_iterator_int it = tab_list[i]->begin(); it != tab_list[i]->end(); ++it) {
         r[counter]= (*it);
         v[counter++]=0;
      }
   }
   for (int i = 0; i<Ng; ++i) delete(tab_list[i]);
   delete[](tab_list);
};

template <typename T>
void inline MaxFlow<T>::save_capacities() {
   for (int i = 0; i<_nzmax; ++i) _copycapacity[i]=_capacity[i];
}
template <typename T>
void inline MaxFlow<T>::save_flow() {
   _copyflow=new T[_nzmax];
   for (int i = 0; i<_nzmax; ++i) _copyflow[i]=_flow[i];
   _copyexcess=new T[_N];
   for (int i = 0; i<_N; ++i) _copyexcess[i]=_excess[i];
}
template <typename T>
void inline MaxFlow<T>::restore_flow() {
   for (int i = 0; i<_nzmax; ++i) _flow[i]=_copyflow[i];
   delete[](_copyflow);
   for (int i = 0; i<_N; ++i) _excess[i]=_copyexcess[i];
   delete[](_copyexcess);
}


template <typename T>
void inline MaxFlow<T>::restore_capacities() {
   for (int i = 0; i<_nzmax; ++i) _capacity[i]=_copycapacity[i];
}

template <typename T>
bool inline MaxFlow<T>::check_flow() {
   list_int tmp;
   for (int i = 0; i<_N; ++i) _seen[i]=false;
   tmp.push_back(_s);
   _seen[_s]=true;

   bool correct=true;
   T total_excess=0;
   T total_excess2=0;
   while (!tmp.empty()) {
      const int node = tmp.front();
      const int ind = _pr_node[node];
      tmp.pop_front();
      if (_excess[node] < 0) {
         cerr << "negative excess: " <<_excess[node]  << " on node " << node << endl;
         correct=false;
      }
      T totalflow=0;
      for (int i = 0; i<_num_edges[node]; ++i) {
         totalflow+=_flow[ind+i];
         if (_flow[ind+i] > _capacity[ind+i]) {
            correct=false;
            cerr << "exceed capacity on node " << node << " to " << _children[ind+i] << ". Flow: " << _flow[ind+i] << ", capacity: " << _capacity[ind+i] << endl;
            total_excess += _flow[ind+i]-_capacity[ind+i];
         }
         if (!_seen[_children[ind+i]]) {
            tmp.push_back(_children[ind+i]);
            _seen[_children[ind+i]]=true;
         }
      }
      if (node != _s && node != _t && abs(totalflow+_excess[node]) > EPSILON_MAXFLOW) {
         cerr << "prb on node " << node << ", excess: " << _excess[node] << ", totalflow: " << totalflow << endl;
      }
      if (node != _s && node != _t)
         total_excess2+=abs(totalflow+_excess[node]);
   }
   return correct;
}

template <typename T>
void inline MaxFlow<T>::reset_flow() {
   memset(_excess,0,_N*sizeof(T));
   memset(_flow,0,_nzmax*sizeof(T));
   _excess[_s]=INFINITY;
}

template <typename T>
void inline MaxFlow<T>::scale_flow(const T scal) {
   for (int i = 0; i<_N; ++i) _excess[i]*=scal;
   for (int i = 0; i<_nzmax; ++i) _flow[i]*=scal;
   _excess[_s]=INFINITY;
}

template <typename T>      
void inline MaxFlow<T>::set_weights(const T lambda) {
   for (int j = 0; j<_num_edges[_s]; ++j) {
      _capacity[_pr_node[_s]+j]=lambda;
   }
};

template <typename T>      
void inline MaxFlow<T>::set_weights(const T* weights, 
      const T lambda) {
   for (int j = 0; j<_num_edges[_s]; ++j) {
      _capacity[_pr_node[_s]+j]=lambda*weights[j];
   }
};

template <typename T>
void inline MaxFlow<T>::print_sink() {
   cerr << "Flow: ";
   for (int j = 0; j<_num_edges[_t]; ++j) {
      cerr << _flow[_reverse_address[_pr_node[_t]+j]] << " ";
   }
   cerr << endl;
   cerr << "Capacity: ";
   for (int j = 0; j<_num_edges[_t]; ++j) {
      cerr << _capacity[_reverse_address[_pr_node[_t]+j]] << " ";
   }
   cerr << endl;
};

template <typename T>
void inline MaxFlow<T>::print_component(const list_int& component) {
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      cerr << "Node: " << *it << endl;
      cerr << "Children: ";
      for (int j = 0; j<_num_edges[*it]; ++j) {
         cerr << _children[_pr_node[*it]+j] << " ";
      }
      cerr << endl;
      cerr << "Flow: ";
      for (int j = 0; j<_num_edges[*it]; ++j) {
         cerr << _flow[_pr_node[*it]+j] << " ";
      }
      cerr << endl;
      cerr << "Capacity: ";
      for (int j = 0; j<_num_edges[*it]; ++j) {
         cerr << _capacity[_pr_node[*it]+j] << " ";
      }
      cerr << endl;
   }
};

template <typename T>
void inline MaxFlow<T>::sub_gradient_aux(const Vector<T>& input, Vector<T>& output, const Vector<T>& weights,
      const int node, list_int& variables, const int Ng) {
   _seen[node]=true;
   const int ind = _pr_node[node];
   const int* children = _children+ind;
   const T* capacity = _capacity+ind;
   for (int i = 0; i<_num_edges[node]; ++i) {
      const int child=children[i];
      if (child != _s && child != _t) {
         if (child < Ng) {
            if (capacity[i] > 0 && !_seen[child]) {
               list_int new_var;
               this->sub_gradient_aux(input,output,weights,child,new_var,Ng);
               variables.fusion(new_var);
            }
         } else {
            variables.push_back(child);
         }
      }  
   }
   T max_abs = 0;
   for (const_iterator_int it = variables.begin(); it != variables.end(); ++it) {
      if (abs<T>(input[*it-Ng]) > max_abs) max_abs=abs<T>(input[*it-Ng]);
   }
   if (max_abs < 1e-15)
      return;
   list_int var_max;
   for (const_iterator_int it = variables.begin(); it != variables.end(); ++it) {
      if (abs<T>(abs<T>(input[*it-Ng])-max_abs) < 1e-15)
         var_max.push_back(*it-Ng);
   }
   T scal = weights[node]/var_max.size();
   for (const_iterator_int it = var_max.begin(); it != var_max.end(); ++it) {
      output[*it] += input[*it] > 0 ? scal : -scal;
   }
};

template <typename T>
void inline MaxFlow<T>::sub_gradient(const Vector<T>& input, Vector<T>& output, const Vector<T>& weights, const int Ng) {
   output.setZeros();
   list_int tmp;
   for (int i = 0; i<_N; ++i) {
      _seen[i]=false;
      if (i < Ng) tmp.push_back(i);
   }
   while (!tmp.empty()) {
      const int node=tmp.front();
      if (!_seen[node]) {
         list_int variables;
         this->sub_gradient_aux(input,output,weights,node,variables,Ng);
      }
      tmp.pop_front();
   }
};

template <typename T>
T inline MaxFlow<T>::norm(const T* variables, T* work, const T* weights, const int Ng, const bool linf) {
   list_int tmp;
   for (int i = 0; i<_N; ++i) {
      _seen[i]=false;
      work[i]=0;
      if (i < Ng) tmp.push_back(i);
   }

   while (!tmp.empty()) {
      const int node = tmp.front();
      if (_seen[node]) {
         tmp.pop_front();
      } else {
         if (node >= Ng && node != _s && node != _t) {
            work[node]= linf ? abs(variables[node-Ng]) : variables[node-Ng]*variables[node-Ng];
            _seen[node]=true;
            tmp.pop_front();
         } else {
            const int ind = _pr_node[node];
            const int* children = _children+ind;
            const T* capacity = _capacity+ind;
            bool all_child=true;
            for (int i = 0; i<_num_edges[node]; ++i) {
               const int child=children[i];
               if (child != _s && child != _t && capacity[i] > 0 && !_seen[child]) {
                  all_child=false;
                  tmp.push_front(children[i]);
               }
            }
            if (all_child) {
               T max_var=0;
               for (int i = 0; i<_num_edges[node]; ++i) {
                  const int child=children[i];
                  if (child != _s && child != _t && capacity[i] > 0) {
                     max_var = linf ? MAX(max_var,work[child]) : max_var + work[child];
                  }
               }
               work[node]=max_var;
               _seen[node]=true;
            }
         }
      }
   }
   T sum=0;
   if (linf) {
      for (int i = 0; i<Ng; ++i) {
         sum+=weights[i]*work[i];
      }
   } else {
      for (int i = 0; i<Ng; ++i) {
         sum+=weights[i]*sqrt(work[i]);
      }
   }
   return sum;
};

template <typename T>
void inline MaxFlow<T>::print_component2(const list_int& component) {
   cerr << "*********Print component***********" << endl;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      cerr << *it <<  " ";
   }
   cerr << endl;
   cerr << "Excess" << endl;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      cerr << _excess[*it] <<  " ";
   }
   cerr << "  " << _excess[_s] << " " << _excess[_t];
   cerr << endl;
   cerr << "Labels" << endl;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      cerr << _labels[*it] <<  " ";
   }
   cerr << "  " << _labels[_s] << " " << _labels[_t];
   cerr << endl;


};


template <typename T>
void inline MaxFlow<T>::print_graph_aux(const int i) {
   if (_seen[i]) return;
   cerr << endl;
   cerr << "Node: " << i << endl;
   _seen[i]=true;
   if (i == _t) return;
   cerr << "Children: ";
   for (int j = 0; j<_num_edges[i]; ++j) {
      cerr << _children[_pr_node[i]+j] << " ";
   }
   cerr << endl;
   cerr << "Capacity: ";
   for (int j = 0; j<_num_edges[i]; ++j) {
      cerr << _capacity[_pr_node[i]+j] << " ";
   }
   cerr << endl;
   cerr << "Flow: ";
   for (int j = 0; j<_num_edges[i]; ++j) {
      cerr << _flow[_pr_node[i]+j] << " ";
   }
   cerr << endl;
   cerr << "Rverse Flow: ";
   for (int j = 0; j<_num_edges[i]; ++j) {
      cerr << _flow[_reverse_address[_pr_node[i]+j]] << " ";
   }
   cerr << endl;

   for (int j = 0; j<_num_edges[i]; ++j) {
      this->print_graph_aux(_children[_pr_node[i]+j]);
   }
};

template <typename T>
inline void MaxFlow<T>::init_split_variables_aux(const int i,int& current,
      Vector<int>& count, list_int** splitted_w, const int Ng, const int Nv) {
   if (_seen[i]  || (i >= Ng && i != _s)) return;
   _seen[i]=true;
   const int ind = _pr_node[i];
   const int* children = _children+ind;
   const T* capacity = _capacity+ind;

   for (int j = 0; j<_num_edges[i]; ++j) {
      if (capacity[j] > 0) {
         this->init_split_variables_aux(children[j],current,count,splitted_w,Ng,Nv);
      }
   }
   if (i != _s) {
      Vector<T> tmp(Nv);
      tmp.setZeros();
      /// rempli colonne current de splitted_w avec les enfants + propres variables
      for (int j = 0; j<_num_edges[i]; ++j) {
         const int child=children[j];
         if (child != _s && child != _t && capacity[j] > 0) {
            if (child >= Ng) {
               tmp[child-Ng]=1.0;
            } else {
               for (const_iterator_int it = splitted_w[count[child]]->begin();
                     it != splitted_w[count[child]]->end(); ++it)
                  tmp[*it]++;
            }
         }
      }
      for (int j = 0; j<tmp.n(); ++j) {
         if (tmp[j]) splitted_w[current]->push_back(j);
      }
      count[i]=current;
      ++current;
   }
};

template <typename T>
void inline MaxFlow<T>::global_relabelling() {
   // global relabelling by reverse breadth first search
   list_int nodes;
   for (int i = 0; i<_N; ++i) _labels[i]=_N;
   for (int i = 0; i<_N; ++i) _seen[i]=false;
   nodes.push_back(_t);
   _seen[_t]=true;
   _labels[_t]=0;
   while (!nodes.empty()) {
      const int node=nodes.front();
      const int* children=_children+_pr_node[node];
      const int* reverse=_reverse_address+_pr_node[node];
      for (int i = 0; i<_num_edges[node]; ++i) {
         const int child=children[i];
         if (!_seen[child] && _capacity[reverse[i]] > _flow[reverse[i]]) {
            _seen[child]=true;
            const int new_label=_labels[node]+1;
            if (new_label != _labels[child] && _excess[child] > EPSILON_MAXFLOW) {
               _active_nodes[new_label]->push_back(child);
               _active[child]=true;
               if (new_label > _current_max_label) 
                  _current_max_label=new_label;
            }
            _labels[child]=new_label;
            nodes.push_back(child);
         }
      }
      nodes.pop_front();
      _active[node]=false;
   }
};

template <typename T>      
void inline MaxFlow<T>::extractConnexComponents(
      std::list< list_int* >& connex_components) {
   /// extract all the connex components for the initialization
   for (int i = 0; i<_N; ++i) _seen[i]=false;
   _seen[_s]=true;
   _seen[_t]=true;
   list_int tmp;
   for (int i = 0; i<_N; ++i) {
      if (!_seen[i]) {
         // create a component
         list_int* component = new list_int();
         tmp.push_back(i);
         while (!tmp.empty()) {
            int node=tmp.front();
            _seen[node]=true;
            component->push_back(node);
            tmp.pop_front();
            const int* children=_children+_pr_node[node];
            for (int i = 0; i<_num_edges[node]; ++i) {
               const int child=children[i];
               if (!_seen[child]) {
                  _seen[child]=true;
                  tmp.push_back(child);
               }
            }
         }
         connex_components.push_back(component);
      }
   }
};

template <typename T>
T inline MaxFlow<T>::project_weighted(const list_int& component, 
      const T* variables_in,T* variables_out,T* work, const T* weights,
      const int Ng) {
   T lambda=0;
   int num_var=0;
   Vector<T> ww(component.size());
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      if (*it < Ng) {
         lambda+=_capacity[_reverse_address[_pr_node[*it]]];
      } else {
         ww[num_var]=T(1.0)/weights[*it-Ng];
         work[num_var++]=variables_in[*it-Ng];
      }
   }
   ww.setn(num_var);
   Vector<T> out;
   Vector<T> in(work,num_var);
   in.l1project_weighted(out,ww,lambda,true);
   T max_flow=0;
   int count=0;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      const int ind = _pr_node[*it];
      if (*it >= Ng) {
         const int nv=*it-Ng;
         variables_out[nv]=out[count];
         const T diff = (variables_in[nv]-variables_out[nv])*ww[count++];
         max_flow+=diff;
         _capacity[ind]=diff;
         if (_flow[ind] > diff) {
            _excess[*it]+=_flow[ind]-diff;
            _flow[ind]=diff;
            _flow[_reverse_address[ind]]=-diff;
         }
         _labels[*it]=1;
      }
   }
   return max_flow;

};

template <typename T>
T inline MaxFlow<T>::project(const list_int& component, 
      const T* variables_in,T* variables_out, T* work,
      const int Ng) {
   /// project on the component, project, update the capacity,
   /// update the preflow,  update variables_out,
   /// update labels 
   /// return the maximum value of the potential flow
   T lambda=0;
   int num_var=0;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      if (*it < Ng) {
         lambda+=_capacity[_reverse_address[_pr_node[*it]]];
      } else {
         work[num_var++]=variables_in[*it-Ng];
      }
   }
   //  PRINT_F(lambda)
   T max_flow=0;
   const T thrs=this->compute_thrs_project_l1(work,num_var,lambda);
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      const int ind = _pr_node[*it];
      if (*it >= Ng) {
         const int nv=*it-Ng;
         variables_out[nv]=MIN(variables_in[nv],thrs);
         const T diff = variables_in[nv]-variables_out[nv];
         max_flow+=diff;
         _capacity[ind]=diff;
         if (_flow[ind] > diff) {
            _excess[*it]+=_flow[ind]-diff;
            _flow[ind]=diff;
            _flow[_reverse_address[ind]]=-diff;
         }
         _labels[*it]=1;
      }
   }
   return max_flow;
};

template <typename T>
T inline MaxFlow<T>::project_box(const list_int& component, 
      const T* variables_in,T* variables_out,T* work, bool& fusion, 
      const int Ng) {
   list_int nodes;
   list_int variables;
   _seen[_s]=true;
   _active[_s]=false;
   T lambda=0;
   int num_nodes=0;
   for (const_iterator_int it=component.begin();
         it !=component.end(); ++it) {
      const int node = *it;
      const int ind = _pr_node[node];
      _active[node]=true;
      _seen[node]=false;
      if (node < Ng) {
         work[node]=_capacity[_reverse_address[ind]];
         nodes.push_back(node);
         _all_nodes[node]=1;
         ++num_nodes;
         lambda+=work[node];
      } else {
         work[node]=0;
         variables.push_back(*it);
      }
   }
   //  PRINT_F(lambda)
   fusion = num_nodes > 1;
   while (!nodes.empty()) {
      const int node = nodes.front();
      if (!_seen[node]) {
         const int ind = _pr_node[node];
         const int* childrens=_children+ind;
         const int* reverse=_reverse_address+ind;
         for (int& i = _all_nodes[node]; i<_num_edges[node]; ++i) {
            const int child = childrens[i];
            if (_active[child] && !_seen[child] 
                  && _capacity[reverse[i]] > 0) {
               nodes.push_front(child);
               break;
            }
         }
         if (_all_nodes[node]==_num_edges[node]) {
            _seen[node]=true;
            for (int i = 1; i<_num_edges[node]; ++i) {
               const int child = childrens[i];
               if (_active[child] && _capacity[ind+i] > 0) {
                  work[child]+=work[node];
               }
            }
            nodes.pop_front();
         }
      } else {
         nodes.pop_front();
      }
   }

   T thrs = INFINITY;
   if (lambda > 0) {
      std::list<T> var;
      for (const_iterator_int it=variables.begin();
            it !=variables.end(); ++it) {
         if (variables_in[*it-Ng] > 0) {
            var.push_back(variables_in[*it-Ng]);
            T diff = variables_in[*it-Ng]-work[*it];
            if (diff > 0)
               var.push_back(-diff);
         }
      }
      var.sort(compare_abs<T>);
      int num=0;
      T sum=0;
      bool br=false;
      T pivot=0;
      for (typename std::list<T>::const_iterator it=var.begin();
            it != var.end(); ++it) {
         pivot=*it;
         sum+= pivot;
         num+= pivot > 0 ? 1 : -1;
         if (sum-abs<T>(*it)*num > lambda) {
            sum-=*it;
            num-= (*it > 0 ? 1 : -1);
            br=true;
            thrs= num==0 ? pivot : (sum-lambda)/num;
            break;
         }
      }
      if (!br) thrs=MAX(0,num==0 ? pivot : (sum-lambda)/num);
   }
   for (const_iterator_int it=variables.begin();
         it !=variables.end(); ++it) {
      variables_out[*it-Ng] = MIN(MAX(thrs,variables_in[*it-Ng]-work[*it]),
            variables_in[*it-Ng]);
   }
   T maxflow=0;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      const int ind = _pr_node[*it];
      if (*it >= Ng) {
         const int nv=*it-Ng;
         const T diff = variables_in[nv]-variables_out[nv];
         maxflow+=diff;
         _capacity[ind]=diff;
         if (_flow[ind] > diff) {
            _excess[*it]+=_flow[ind]-diff;
            _flow[ind]=diff;
            _flow[_reverse_address[ind]]=-diff;
         }
         _labels[*it]=1;
      }
   }
   return maxflow;
};


template <typename T>
T inline MaxFlow<T>::flow_component(const list_int& component, const int Ng) {
   /// do relabelling and update list of active nodes
   T max_flow=0;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      if (*it >= Ng) {
         max_flow+=_flow[_pr_node[*it]];
      }
   }
   /// return the amount of flow
   return max_flow;
};

/*template <typename T>
  bool inline MaxFlow<T>::splitComponent2(const list_int& component,
  std::list< list_int* >& connex_components,const int Ng, bool* positive, const bool addpos) {

  }*/

template <typename T>
bool inline MaxFlow<T>::splitComponent(const list_int& component,
      std::list< list_int* >& connex_components,const int Ng, bool* positive, const bool addpos) {
   /// cut the component into connex components, and add them to connex_components
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      _seen[*it]=false;
      positive[*it]=false;
   }
   int num_comp=0;
   _seen[_s]=true;
   _seen[_t]=true;
   positive[_s]=true;
   positive[_t]=true;
   list_int tmp;
   /// make the "positive part of the graph"
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      if (!positive[*it] && _excess[*it] > EPSILON_MAXFLOW) {
         /// start new component, track from where the excess can come
         tmp.push_back(*it);
         positive[*it]=true;
         while (!tmp.empty()) {
            int node=tmp.front();
            tmp.pop_front();
            const int ind=_pr_node[node];
            const int* children=_children+ind;
            const T* flow=_flow+ind;
            const T* capacity=_capacity+ind;
            for (int i = 0; i<_num_edges[node]; ++i) {
               const int child=children[i];
               if (!_seen[child] && !positive[child] && (flow[i] < capacity[i]-EPSILON_MAXFLOW)) {
                  positive[child]=true;
                  tmp.push_back(child);
               }
            }
         }
      }
   }
   /// update from the source
   /*tmp.push_back(_s);
     while (!tmp.empty()) {
     int node=tmp.front();
     tmp.pop_front();
     const int ind=_pr_node[node];
     const int* children=_children+ind;
     const T* flow=_flow+ind;
     const T* capacity=_capacity+ind;
     for (int i = 0; i<_num_edges[node]; ++i) {
     const int child=children[i];
     if (!_seen[child] && !positive[child] && (flow[i] < capacity[i]-EPSILON_MAXFLOW)) {
     positive[child]=true;
     tmp.push_back(child);
     }
     }
     }*/

   /// extract the connex components of the positive part
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      if (positive[*it] && !_seen[*it]) {
         list_int* new_component = new list_int();
         /// start new component, track from where the excess can come
         tmp.push_back(*it);
         _seen[*it]=true;
         while (!tmp.empty()) {
            int node=tmp.front();
            new_component->push_back(node);
            tmp.pop_front();
            const int ind=_pr_node[node];
            const int* children=_children+ind;
            for (int i = 0; i<_num_edges[node]; ++i) {
               const int child=children[i];
               if (!positive[child] && child != _t) {
                  _capacity[ind+i]=_capacity[ind+i] > 0 ? -0.5 : 0;
               }
               if (positive[child] && !_seen[child]) {
                  _seen[child]=true;
                  tmp.push_back(child);
               }
            }
         }
         if (addpos) {
            connex_components.push_back(new_component);
         } else {
            delete(new_component);
         }
         ++num_comp;
      }
   }
   /// extract the connex components of the negative part
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      if (!positive[*it] && !_seen[*it]) {
         list_int* new_component = new list_int();
         /// start new component, track from where the excess can come
         tmp.push_back(*it);
         _seen[*it]=true;
         while (!tmp.empty()) {
            int node=tmp.front();
            new_component->push_back(node);
            tmp.pop_front();
            const int ind=_pr_node[node];
            const int* children=_children+ind;
            for (int i = 0; i<_num_edges[node]; ++i) {
               const int child=children[i];
               if (positive[child] && child != _t) {
                  //_capacity[ind+i]=-0.5;
                  _capacity[ind+i]=_capacity[ind+i] > 0 ? -0.5 : 0;
               }
               if (!positive[child] && !_seen[child]) {
                  _seen[child]=true;
                  tmp.push_back(child);
               }
            }
         }
         connex_components.push_back(new_component);
         ++num_comp;
      }
   }
   if (num_comp == 1 && connex_components.size() != 0) {
      list_int* comp = connex_components.back();
      delete(comp);
      connex_components.pop_back();
   }
   return num_comp > 1;
   // cout << "Number of new component: " << num_comp << endl;
};

template <typename T>
void inline MaxFlow<T>::reset_component(const list_int& component) {
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      const int ind = _pr_node[*it];
      _excess[*it]=0;
      for (int i = 0; i<_num_edges[*it]; ++i) {
         _flow[i+ind]=0;
         _flow[_reverse_address[i+ind]]=0;
      }
   }
};

template <typename T>
void inline MaxFlow<T>::perform_maxflow_component(const list_int& component) {
   tglobal3.start();
   int size_component=component.size();
   const int max_label=size_component+2;
   /// discharge the source and relabel
   this->component_relabelling(component,max_label,true);
#ifdef VERBB
   PRINT_I(_current_max_label)
      this->print_component2(component);
   this->print_component(component);
   stop();
#endif
   /// perform max flow
   int counter=1;


   while (_current_max_label > 0 || !_active_nodes[0]->empty()) { 
#ifdef VERBB
      this->print_component2(component);
      stop();
#endif
      if (global_heuristic && (counter % (size_component+1)) == 0) {
         this->component_relabelling(component,max_label,false);
         ++counter;
      } else {
         if (_active_nodes[_current_max_label]->empty()) {
            _current_max_label--;
#ifdef VERBB
            cerr << "current max label decreased to " << _current_max_label << endl;
#endif
         } else {
            const int current_node=_active_nodes[_current_max_label]->front();
            _active_nodes[_current_max_label]->pop_front();
            _active[current_node]=false;
            if (_excess[current_node] > EPSILON_MAXFLOW) {
               this->discharge(component,current_node,max_label);
               if (_excess[current_node] > EPSILON_MAXFLOW && _labels[current_node] < max_label) {
                  _active_nodes[_labels[current_node]]->push_back(current_node);
                  _active[current_node]=true;
                  if (_labels[current_node] > _current_max_label) {
                     _current_max_label=_labels[current_node];
                  }
               }
            } else {
               _excess[current_node]=0;
            }
            ++counter;
         }
      }
   }
#ifdef VERBB
   cerr << "END max flow" << endl;
   this->print_excess();
   stop();
#endif
   tglobal3.stop();
};

template <typename T>
void inline MaxFlow<T>::gap_relabelling(const list_int& component, const int gap, const int max_label) {
#ifdef VERBB
   cerr << "Gap relabelling " << gap << endl;
#endif
   if (tglobal2.getElapsed() > 0.1*tglobal3.getElapsed()) return;
   tglobal2.start();
   num_gap_relabels++;
   for (const_iterator_int it = component.begin(); it != component.end(); ++it) {
      if (_labels[*it] > gap) {
         _labels[*it]=max_label;
      }
   }
   for (int i = gap; i<max_label; ++i) {
      _all_nodes[i]=0;
   }
   tglobal2.stop();
};

template <typename T>
void inline MaxFlow<T>::component_relabelling(const list_int& component,
      const int max_label, const bool force) {
   tglobal1.start();
   if (!force && (tglobal1.getElapsed() > 0.1*tglobal3.getElapsed())) return;
   for (int i = 0; i<=component.size(); ++i)
      _active_nodes[i]->clear();
   if (gap_heuristic)
      for (int i = 0; i<=component.size(); ++i)
         _all_nodes[i]=0;
   _current_max_label=0;
   num_global_relabels++;
   /// relabel component, with warm restart
   list_int nodes;
   _labels[_t]=0;
   _all_nodes[0]++;
   _labels[_s]=max_label;
   _seen[_t]=true;
   _active[_t]=false;
   _seen[_s]=true;
   _active[_s]=false;
   //  _seen[_s]=false;
   //  _active[_s]=true;

   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      const int ind = _pr_node[*it];
      const int first_child=_children[ind];
      if (first_child==_t && _flow[ind] < _capacity[ind]) {
         _labels[*it]=1;
         nodes.push_back(*it);
         if (_excess[*it] > EPSILON_MAXFLOW) {
            _active_nodes[1]->push_back(*it);
            _current_max_label=1;
            _active[*it]=true;
         } else {
            _active[*it]=false;
         }
         if (gap_heuristic)
            _all_nodes[1]++;
         _seen[*it]=true;
      } else {
         /// discharge source
         if (first_child == _s && force) {
            const T delta = _capacity[_reverse_address[ind]] - _flow[_reverse_address[ind]];
            if (delta > 0) {
               _excess[*it] += delta;
               _flow[_reverse_address[ind]]=_capacity[_reverse_address[ind]];
            }
         }
         _seen[*it]=false;
         _active[*it]=false;
         _labels[*it]=max_label;
      }
   }
   while (!nodes.empty()) {
      const int node=nodes.front();
      const int* children=_children+_pr_node[node];
      const int* reverse=_reverse_address+_pr_node[node];
      for (int i = 0; i<_num_edges[node]; ++i) {
         const int child=children[i];
         if (!_seen[child] && _capacity[reverse[i]] > _flow[reverse[i]]) {
            _seen[child]=true;
            const int new_label=_labels[node]+1;
            if (new_label != _labels[child] && _excess[child] > EPSILON_MAXFLOW) {
               _active_nodes[new_label]->push_back(child);
               _active[child]=true;
               if (new_label > _current_max_label) 
                  _current_max_label=new_label;
            }
            _labels[child]=new_label;
            if (gap_heuristic)
               _all_nodes[new_label]++;
            nodes.push_back(child);
         }
      }
      nodes.pop_front();
   }
   tglobal1.stop();
   /*   this->print_graph();
        this->print_excess();
        this->print_labels();
        stop();*/
};

template <typename T>
void inline MaxFlow<T>::update_capacities(const list_int& component, T* work) {
   list_int comp_nodes;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      const int ind = _pr_node[*it];
      const int first_child=_children[ind];
      _all_nodes[*it]=0;
      _active[*it]=true;
      if (first_child == _t) {
         _seen[*it]=true;
         work[*it]=_capacity[ind];
      } else {
         _seen[*it]=false;
         comp_nodes.push_back(*it);
      }
   }
   list_int nodes;
   while (!comp_nodes.empty()) {
      const int new_node=comp_nodes.front();
      comp_nodes.pop_front();
      if (!_seen[new_node]) {
         nodes.push_back(new_node);
         while (!nodes.empty()) {
            const int node = nodes.front();
            _seen[node]=true;
            const int ind = _pr_node[node];
            const int* children=_children+ind;
            for ( ; _all_nodes[node] < _num_edges[node]; ++_all_nodes[node]) {
               const int child=children[_all_nodes[node]];
               if (_active[child] && !_seen[child] &&_capacity[ind+_all_nodes[node]] > 0) {
                  nodes.push_front(child);
                  break;
               }
            }
            if (_all_nodes[node]==_num_edges[node]) {
               work[node]=0;
               for (int i = 0; i < _num_edges[node]; ++i) {
                  const int child=children[i];
                  if (_active[child] && _capacity[ind+i] > 0) {
                     if (work[child] > 0) {
                        work[node]+=work[child];
                        _capacity[ind+i] = MAX(_flow[ind+i],work[child]);
                     } else {
                        _capacity[ind+i]=-2;
                     }
                  }
               }
               nodes.pop_front();
            } 
         }
      }
   }
}
template <typename T>
void inline MaxFlow<T>::set_capacities_variables(const T* cap, const int Nv,
      const int Ng) {
   for (int i = 0; i<Nv; ++i) {
      const int ind = _pr_node[Ng+i];
      _capacity[ind]=abs(cap[i]);
   }
};

template <typename T>
void inline MaxFlow<T>::set_capacities_groups(const list_int& component,
      const Vector<T>& weights,const T lambda, const int Ng) {
   for (const_iterator_int it = component.begin(); it != component.end();
         ++it) {
      if (*it < Ng) {
         _capacity[_reverse_address[_pr_node[*it]]]=lambda*weights[*it];
      }
   }
};


template <typename T>
void inline MaxFlow<T>::restore_capacities(const list_int& component) {
   /// relabel component, with warm restart
   list_int nodes;
   _seen[_t]=true;
   _seen[_s]=true;
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      _seen[*it]=false;
   }
   for (const_iterator_int it=component.begin();
         it != component.end(); ++it) {
      const int* children=_children+_pr_node[*it];
      T* capacity=_capacity+_pr_node[*it];
      for (int i = 0; i<_num_edges[*it]; ++i) {
         const int child=children[i];
         if (!_seen[child] && (capacity[i] > 0 || capacity[i] < -1))
            capacity[i]=INFINITY;
      }
   }
};




template <typename T>
T inline MaxFlow<T>::compute_thrs_project_l1(T* X, const int n, const T lambda) {
   if (lambda==0) return INFINITY;
   T* prU = X;
   T sum = 0;
   int sum_card = n;
   for (int i = 0; i<sum_card; ++i) {
      if (X[i]) {
         sum += X[i];
      } else {
         swap(X[i],X[--sum_card]);
         --i;
      }
   }
   if (sum < lambda) {
      memset(X,0,sum_card*sizeof(T));
      return 0;
   }
   int sizeU = sum_card;
   sum_card = 0;
   sum=0;

   while (sizeU > 0) {
      // put the pivot in prU[0]
      swap(prU[0],prU[sizeU/2]);
      int sizeG=1;
      T sumG=prU[0];

      for (int i = 1; i<sizeU; ++i) {
         if (prU[i] >= prU[0]) {
            sumG += prU[i];
            swap(prU[sizeG++],prU[i]);
         }
      }

      T new_sum=sum+sumG;
      int new_card=sum_card+sizeG;
      if (new_sum - prU[0]*new_card <= lambda) {
         sum_card = new_card;
         sum = new_sum;
         prU +=sizeG;
         sizeU -= sizeG;
      } else {
         ++prU;
         sizeU = sizeG-1;
      }
   }
   return MAX(0,(sum-lambda)/sum_card);
};

template <typename T> class Graph {
   public:
      Graph();
      ~Graph();

      void inline create_graph(const int Nv, const int Ng,
            T* weights, mwSize* var_ir, mwSize* var_jc);
      void inline create_graph(const int Nv, const int Ng,
            T* weights, mwSize* gv_ir, mwSize* gv_jc, mwSize* gg_ir, mwSize* gg_jc);

      void inline proximal_operator(const T* variables_in, T* variables_out,const bool clever = false, const T* weights = NULL);
      void inline save_capacities() { _maxflow->save_capacities(); };
      void inline restore_capacities() { _maxflow->restore_capacities(); };
      void inline save_flow() { _maxflow->save_flow(); };
      void inline restore_flow() { _maxflow->restore_flow(); };
      void inline reset_flow() { _maxflow->reset_flow(); };
      void inline scale_flow(const T scal) { _maxflow->scale_flow(scal); };
      void inline set_weights(const T lambda) {
         _maxflow->set_weights(lambda); };
      void inline set_weights(const T* weights,
            const T lambda) {
         _maxflow->set_weights(weights,lambda); };
      void inline print() { _maxflow->print_graph(); };
      T inline norm(const T* variables, T* work, const T* weights, const bool linf=true) { return _maxflow->norm(variables,work,weights,_Ng,linf); };
      T inline dual_norm_inf(const Vector<T>& input, const Vector<T>& weights);
      void inline sub_gradient(const Vector<T>& input, Vector<T>& output, const Vector<T>& weights) {
         _maxflow->sub_gradient(input,output,weights,_Ng);
      }
      inline void init_split_variables(SpMatrix<T>& splitted_w) {
         _maxflow->init_split_variables(splitted_w,_Ng,_Nv);
      };


   private:
      int _Nv;
      int _Ng;
      T* _weights; // size Ng
      MaxFlow<T>* _maxflow;
};

template <typename T>
Graph<T>::Graph() {
   _Nv=0;
   _Ng=0;
   _weights=NULL;
   _maxflow=NULL;
};

template <typename T>
Graph<T>::~Graph() {
   delete[](_weights);
   delete(_maxflow);
};

template <typename T>
void inline Graph<T>::create_graph(const int Nv, const int Ng,
      T* weights, mwSize* var_ir, mwSize* var_jc) {
   _Nv=Nv;
   _Ng=Ng;
   _weights=new T[_Ng];
   for (int i = 0; i<_Ng; ++i) _weights[i]=weights[i];
   const int N = _Ng+_Nv+2;
   int* num_edges=new int[N];
   for (int i = 0; i<N; ++i) num_edges[i]=1;
   for (int i = 0; i<Ng; ++i) {
      for (int j = var_jc[i]; j<var_jc[i+1]; ++j) {
         num_edges[i]++;
         num_edges[Ng+var_ir[j]]++;
      }
   }
   const int s=_Ng+_Nv;
   const int t=_Ng+_Nv+1;
   num_edges[s]=_Ng;
   num_edges[t]=_Nv;
   _maxflow=new MaxFlow<T>(N, num_edges, s, t);

   for (int i = 0; i<_Ng; ++i)
      _maxflow->add_edge(s,i,_weights[i],0);
   for (int i = 0; i<_Nv; ++i)
      _maxflow->add_edge(_Ng+i,t,INFINITY,0);
   for (int i = 0; i<_Ng; ++i) {
      for (int j = var_jc[i]; j<var_jc[i+1]; ++j) {
         _maxflow->add_edge(i,_Ng+static_cast<int>(var_ir[j]),_weights[i],0);
      }
   }
   _maxflow->save_capacities();
   delete[](num_edges);
};

template <typename T>
void inline Graph<T>::create_graph(const int Nv, const int Ng,
      T* weights, mwSize* gv_ir, mwSize* gv_jc, mwSize* gg_ir, mwSize* gg_jc) {
   _Nv=Nv;
   _Ng=Ng;
   _weights=new T[_Ng];
   for (int i = 0; i<_Ng; ++i) _weights[i]=weights[i];
   const int N = _Ng+_Nv+2;
   int* num_edges=new int[N];
   for (int i = 0; i<N; ++i) num_edges[i]=1;
   for (int i = 0; i<Ng; ++i) {
      for (int j = gv_jc[i]; j<static_cast<int>(gv_jc[i+1]); ++j) {
         num_edges[i]++;
         num_edges[Ng+gv_ir[j]]++;
      }
   }
   for (int i = 0; i<Ng; ++i) {
      for (int j = gg_jc[i]; j<static_cast<int>(gg_jc[i+1]); ++j) {
         if (i != static_cast<int>(gg_ir[j])) {
            num_edges[i]++;
            num_edges[gg_ir[j]]++;
         }
      }
   }

   const int s=_Ng+_Nv;
   const int t=_Ng+_Nv+1;
   num_edges[s]=_Ng;
   num_edges[t]=_Nv;

   _maxflow=new MaxFlow<T>(N, num_edges, s, t);

   for (int i = 0; i<_Ng; ++i)
      _maxflow->add_edge(s,i,_weights[i],0);
   for (int i = 0; i<_Nv; ++i)
      _maxflow->add_edge(_Ng+i,t,0,0);
   for (int i = 0; i<_Ng; ++i) {
      for (int j = gv_jc[i]; j<static_cast<int>(gv_jc[i+1]); ++j) {
         _maxflow->add_edge(i,_Ng+static_cast<int>(gv_ir[j]),INFINITY,0);
      }
   }
   for (int i = 0; i<_Ng; ++i) {
      for (int j = gg_jc[i]; j<static_cast<int>(gg_jc[i+1]); ++j) {
         if (i != static_cast<int>(gg_ir[j])) {
            _maxflow->add_edge(i,static_cast<int>(gg_ir[j]),INFINITY,0);
         }
      }
   }
   _maxflow->save_capacities();
   delete[](num_edges);
};

template <typename T>
T inline Graph<T>::dual_norm_inf(const Vector<T>& input,
      const Vector<T>& weights) {
   //  input.print("input");
   Timer time, time2;
   time.start();
   T* work = new T[_Nv+_Ng+2];
   bool* positive = new bool[_Ng+_Nv+2];
   _maxflow->set_capacities_variables(input.rawX(),_Nv,_Ng);
   std::list< list_int* > connex_components;
   _maxflow->extractConnexComponents(connex_components);
   _maxflow->deactivate();
   T tau = 0;
   int num=0;
   long num1=0;
   long num2=0;
   long num3=0;
   long num4=0;

   while (!connex_components.empty()) {
      ++num;
      list_int* component=connex_components.front();
      connex_components.pop_front();
      if (component->size() != 1) {
         // Compute budget and set up input capacities
         T sum_variables=0;
         T sum_weights=0;
         int size_list=0;
         for (const_iterator_int it = component->begin();
               it != component->end(); ++it) {
            if (*it < _Ng) {
               sum_weights+=weights[*it];
               ++size_list;
            } else {
               sum_variables+=abs<T>(input[*it-_Ng]);
            }
         }
         tau = MAX(tau,sum_variables/sum_weights);
         _maxflow->set_capacities_groups(*component,weights,tau,_Ng);
         if (cap_heuristic)
            _maxflow->update_capacities(*component,work);
         num_relabels=0;
         num_global_relabels=0;
         num_gap_relabels=0;
         _maxflow->perform_maxflow_component(*component);
         num1+=num_relabels;
         num2+=num_global_relabels;
         num3+=component->size();
         num4+=num_gap_relabels;
         T flow=_maxflow->flow_component(*component,_Ng);
         _maxflow->restore_capacities(*component);
         if (flow < (sum_variables-EPSILON_MAXFLOW)) {
            _maxflow->splitComponent(*component,connex_components,_Ng, positive,false);
         }
         _maxflow->deactivate(*component);
      }
      delete(component);
   }
   /*cerr << "num_comp: " << num << endl;
     cerr << "size_component: " << num3 << endl;
     cerr << "num_relabels: " << num1 << endl;
     cerr << "global: " << num2 << endl;
     cerr << "gap:: " << num4 << endl;
     cerr << "Time global" << endl;
     cerr << "Time dual" << endl;
     time.printElapsed();*/

   delete[](work);
   delete[](positive);
   return tau;
};

template <typename T>
void inline Graph<T>::proximal_operator(const T* variables_in, T* variables_out,bool clever, const T* weights) {

   tglobal1.reset();
   tglobal2.reset();
   tglobal3.reset();
   Timer tprox;
   tprox.start();
#ifdef VERB2
   PRINT_I(_Nv)
      PRINT_I(_Ng)
      PRINT_I(_maxflow->nzmax())
#endif
      cap_heuristic = true;
   global_heuristic = true;
   gap_heuristic = true;
   std::list< list_int* > connex_components;
   _maxflow->extractConnexComponents(connex_components);
   T* work = new T[_Nv+_Ng+2];
   T* variables_bis = new T[_Nv];
   for (int i = 0; i<_Nv; ++i) variables_bis[i]=abs<T>(variables_in[i]);
   for (int i = 0; i<_Nv; ++i) variables_out[i]=variables_bis[i];

   /*  cerr << "var out" << endl;
       for (int i = 0; i<_Nv; ++i)
       cerr << variables_out[i] << " ";
       cerr << endl;*/
   bool* positive = new bool[_Ng+_Nv+2];
   _maxflow->deactivate();
   T flow_missed=0;
   int num=1;
   long num1=0;
   long num2=0;
   long num3=0;
   long num4=0;

   Timer tsplit, tproj, tcap;

   while (!connex_components.empty()) {
      list_int* component=connex_components.front();
      connex_components.pop_front();
      if (component->size() != 1) {
         bool fusion=true;
         T budget;

         tproj.start();
         if (weights) {
            budget=_maxflow->project_weighted(*component,variables_bis,variables_out,work,weights,_Ng);
         } else {
            if (clever) {
               budget=_maxflow->project_box(*component,variables_bis,variables_out,work,fusion,_Ng);
            } else {
               budget=_maxflow->project(*component,variables_bis,variables_out,work,_Ng);
            }
         }
         tproj.stop();
         ++num;
         if (budget > EPSILON_MAXFLOW && fusion) {
            /// At this point, the vector _maxflow->_seen is set to true.
            tcap.start();
            if (cap_heuristic)
               _maxflow->update_capacities(*component,work);
            tcap.stop();
            num_relabels=0;
            num_global_relabels=0;
            num_gap_relabels=0;
            _maxflow->perform_maxflow_component(*component);
            num1+=num_relabels;
            num2+=num_global_relabels;
            num3+=component->size();
            num4+=num_gap_relabels;
            _maxflow->restore_capacities(*component);
            T flow=_maxflow->flow_component(*component,_Ng);
            if (abs<T>(budget-flow)/budget > EPSILON_MAXFLOW) {
               /// At this point, the vector _maxflow->_seen is set to true for all nodes not in component
               tsplit.start();
               if (!_maxflow->splitComponent(*component,connex_components,_Ng, positive,true)) {
                  flow_missed+=abs<T>(budget-flow);
               }
               tsplit.stop();
            }
         }
#ifdef VERB3
         if (component->size() > 100000) {
            PRINT_I(component->size())
               cerr << "num_comp: " << num << endl;
            cerr << "size_component: " << num3 << endl;
            cerr << "num_relabels: " << num1 << endl;
            cerr << "global: " << num2 << endl;
            cerr << "gap:: " << num4 << endl;
            tglobal1.printElapsed();
            tglobal2.printElapsed();
            tglobal3.printElapsed();
            tcap.printElapsed();
            tproj.printElapsed();
            tsplit.printElapsed();
         }
#endif
         _maxflow->deactivate(*component);
      } 
      delete(component);
   }
#ifdef VERB2
   cerr << "num_comp: " << num << endl;
   cerr << "size_component: " << num3 << endl;
   cerr << "num_relabels: " << num1 << endl;
   cerr << "global: " << num2 << endl;
   cerr << "gap:: " << num4 << endl;
   cerr << "Time global" << endl;
   tglobal1.printElapsed();
   tglobal2.printElapsed();
   tglobal3.printElapsed();
   tcap.printElapsed();
   tproj.printElapsed();
   tsplit.printElapsed();
   tprox.printElapsed();
#endif

   /*  cerr << "var out" << endl;
       for (int i = 0; i<_Nv; ++i)
       cerr << variables_out[i] << " ";
       cerr << endl;*/
   for (int i = 0; i<_Nv; ++i) variables_out[i] = variables_in[i] >= 0 ? MAX(variables_out[i],0) : -MAX(variables_out[i],0);

   delete[](positive);
   delete[](variables_bis);
   delete[](work);
};

template <typename T> struct GraphStruct { 
   mwSize* gv_ir;
   mwSize* gv_jc;
   mwSize* gg_ir;
   mwSize* gg_jc;
   int Nv;
   int Ng;
   T* weights;
};

template <typename T> struct TreeStruct { 
   int* own_variables;
   int* N_own_variables;
   T* weights;
   mwSize* groups_ir;
   mwSize* groups_jc; 
   int Nv;
   int Ng;
};



#endif
