// Full paired-box benchmark. See the repository README for build and audit commands.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
extern "C" {
struct lprec;
lprec* make_lp(int,int);
void delete_lp(lprec*);
void set_verbose(lprec*,int);
unsigned char set_add_rowmode(lprec*,unsigned char);
unsigned char add_constraint(lprec*,double*,int,double);
unsigned char set_bounds(lprec*,int,double,double);
unsigned char set_obj_fn(lprec*,double*);
void set_maxim(lprec*);
int solve(lprec*);
unsigned char get_variables(lprec*,double*);
}
using V=std::vector<double>;
using Matrix=std::vector<V>;
using Labels=std::vector<int>;
using Clock=std::chrono::steady_clock;
std::string context;
#define require(ok) do {if(!(ok))throw std::runtime_error(context+" failed at "+std::to_string(__LINE__)+" "+#ok);} while(false)
struct LP {
    lprec* ptr;
    explicit LP(int cols):ptr(make_lp(0,cols)){require(ptr!=nullptr);set_verbose(ptr,0);}
    ~LP(){delete_lp(ptr);}
    LP(const LP&)=delete;
};
Matrix fixture(int n,int seed,const std::string& family){
    int r=n/2,d=n-r;Matrix m(r,V(d));std::mt19937 rng(seed);
    for(int i=0;i<r;++i)for(int j=0;j<d;++j){
        int v=int(rng()%3)-1;
        if(family=="mixed")m[i][j]=v;
        else if(family=="split")m[i][j]=(i<r/2)==(j<d/2)?(i<r/2?1:-1)*(v!=0):0;
        else if(family=="chain")m[i][j]=i==j?1:(j==i+1?-2:0);
        else throw std::runtime_error("unknown family");
    }
    return m;
}
// Full paired basis L B_*, with L=(I+11^T/(2N))/2. Use both signed slots.
Labels full_box(const Matrix& m,bool batch){
    int r=m.size(),d=m[0].size(),n=r+d,cols=(batch?4:3)*n;
    Matrix base(n,V(2*n));
    for(int i=0;i<r;++i){base[i][i]=1;for(int j=0;j<d;++j)base[i][r+j]=m[i][j];}
    for(int j=0;j<d;++j){base[r+j][n+r+j]=1;for(int i=0;i<r;++i)base[r+j][n+i]=-m[i][j];}
    LP lp(cols);double tau=1./(32*n);
    for(int j=0;j<n;++j)require(set_bounds(lp.ptr,1+j,-6,6));
    for(int j=0;j<2*n;++j)require(set_bounds(lp.ptr,1+n+j,0,2));
    require(set_add_rowmode(lp.ptr,1));
    for(int j=0;j<2*n;++j){
        double shift=0;for(int k=0;k<n;++k)shift+=base[k][j]/(2*n);
        V row(cols+1);for(int k=0;k<n;++k)row[1+k]=(base[k][j]+shift)/2;
        row[1+n+j]=-1;require(add_constraint(lp.ptr,row.data(),3,0));
    }
    if(batch) for(int j=0;j<n;++j){
        V row(cols+1);row[1+n+j]=-1;row[1+3*n+j]=1;
        require(add_constraint(lp.ptr,row.data(),1,0));require(set_bounds(lp.ptr,1+3*n+j,0,tau));
    }
    require(set_add_rowmode(lp.ptr,0)); Labels ans(n);
    if(batch){
        V obj(cols+1),sol(cols);std::fill(obj.begin()+1+3*n,obj.end(),1);
        require(set_obj_fn(lp.ptr,obj.data()));set_maxim(lp.ptr);
        require(solve(lp.ptr)==0 && get_variables(lp.ptr,sol.data()));
        for(int j=0;j<n;++j)ans[j]=sol[3*n+j]>tau/2;
    }else for(int j=0;j<n;++j){
        require(set_bounds(lp.ptr,1+n+j,.25,2));
        int status=solve(lp.ptr);require(status==0 || status==2);ans[j]=(status==0);
        require(set_bounds(lp.ptr,1+n+j,0,2));
    }
    return ans;
}
int main(int argc,char** argv) try {
    if(argc==1){
        std::cout<<"N,family,seed,rep,coordinate_ms,batch_ms,reference_mask,batch_mask\n";
        for(int n:{8,16,32,64})for(auto family:{"mixed","split","chain"})for(int seed:{11,29,47}){
            Matrix m=fixture(n,seed,family);context=std::to_string(n)+" "+family;
            for(int rep=-1;rep<5;++rep){
                V ms(2);Labels result[2];
                for(int j=0;j<2;++j){
                    int mode=(rep+j+2)%2;auto start=Clock::now();
                    result[mode]=full_box(m,mode);
                    ms[mode]=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
                }
                if(rep>=0){
                    std::cout<<n<<','<<family<<','<<seed<<','<<rep<<','<<ms[0]<<','<<ms[1];
                    for(auto& labels:result){std::cout<<',';for(int l:labels)std::cout<<l;}
                    std::cout<<'\n';
                }
            }
        }
        return 0;
    }

    if(argc==2 && std::string(argv[1])=="--matrices"){
        for(int n:{8,16,32,64})for(auto family:{"mixed","split","chain"})for(int seed:{11,29,47}){
            Matrix m=fixture(n,seed,family);
            std::cout<<n<<' '<<family<<' '<<seed;
            for(const auto& row:m)for(double a:row)std::cout<<' '<<int(a);
            std::cout<<'\n';
        }
        return 0;
    }
    throw std::runtime_error("usage: paired [--matrices]");
} catch(const std::exception& ex){std::cerr<<ex.what()<<"\n";return 1;}
