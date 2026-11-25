"""Independent integer-simplex visibility oracle. Run on paired --matrices output and the benchmark CSV."""
from concurrent.futures import ProcessPoolExecutor
import sys
import csv

class Cube:
    def __init__(self, forms):
        self.forms=forms
        self.r=len(forms[0])
        # All head forms lie in [0,1], generator coefficients are nonnegative.
        # Slack dictionary x_basic=(b-A*x_nonbasic)/den.
        self.rows=[[*a,1] for a in forms]+[[-v for v in a]+[0] for a in forms if min(a)<0]
        self.den=1
        self.nonbasic=list(range(self.r))
        self.basic=list(range(self.r,self.r+len(self.rows)))

    def maximum(self, objective):
        cost=list(objective)+[0]*len(self.basic)
        seen=set()
        bland=False
        while True:
            active=[(row,cost[b]) for row,b in zip(self.rows,self.basic) if cost[b]]
            gain=[self.den*cost[j]-sum(row[k]*c for row,c in active) for k,j in enumerate(self.nonbasic)]
            candidates=[k for k,v in enumerate(gain) if v>0]
            if not candidates:
                return sum(row[-1]*c for row,c in active),self.den
            key=tuple(self.basic)
            if key in seen:
                bland=True
            seen.add(key)
            col=min(candidates,key=lambda k:(0 if bland else -gain[k],self.nonbasic[k]))
            leave=None
            for i,row in enumerate(self.rows):
                if row[col]<=0:
                    continue
                if leave is None:
                    leave=i
                    continue
                other=self.rows[leave]
                diff=row[-1]*other[col]-other[-1]*row[col]
                if diff<0 or (diff==0 and self.basic[i]<self.basic[leave]):
                    leave=i
            assert leave is not None
            pivot=self.rows[leave]
            newden=pivot[col]
            den=self.den
            for i,row in enumerate(self.rows):
                if i==leave:
                    continue
                mul=row[col]
                if mul==0 and newden==den:
                    continue
                # Fraction-free exchange by the determinant identity.
                new=[(a*newden-b*mul)//den for a,b in zip(row,pivot)]
                assert all((a*newden-b*mul)%den==0 for a,b in zip(row,pivot))
                new[col]=-mul
                assert new[-1]>=0
                self.rows[i]=new
            pivot[col]=den
            self.den=newden
            self.nonbasic[col],self.basic[leave]=self.basic[leave],self.nonbasic[col]

def check(line):
    n,fam,seed,*fields=line.split()
    n=int(n);r=n//2
    data=list(map(int,fields[:r*r]))
    m=[data[j:j+r] for j in range(0,r*r,r)]
    eye=[[int(i==j) for j in range(r)] for i in range(r)]
    primal=Cube(eye+list(zip(*m)))
    dual=Cube([[-v for v in row] for row in m]+eye)
    labels=[[int(c) for c in s] for s in fields[r*r:]]
    pos=neg=0
    failures=[0]*len(labels)
    baseline_diff=0
    for i in range(n):
        a,b=primal.maximum(primal.forms[i])
        c,d=dual.maximum(dual.forms[i])
        assert 0<=a<=b and 0<=c<=d and a*c==0
        baseline_diff+=(labels[0][i]!=(8*a>=b))
        required=1 if 3*a>=b else (0 if 3*c>=d else None)
        if required is not None:
            pos+=(required==1)
            neg+=(required==0)
            for k in range(len(labels)):
                failures[k]+=(labels[k][i]!=required)
    return (n,fam,seed,pos,neg,baseline_diff,*failures)

if __name__=='__main__':
    with open(sys.argv[1]) as f:
        inputs=f.readlines()
    paired={}
    with open(sys.argv[2]) as f:
        for row in csv.DictReader(f):
            key=(row['N'],row['family'],row['seed'])
            masks=(row['reference_mask'],row['batch_mask'])
            assert masks == paired.setdefault(key,masks)
    assert len(inputs)==len(paired)==36
    records=[]
    for line in inputs:
        key=tuple(line.split()[:3])
        masks=paired.pop(key)
        assert all(len(mask)==int(key[0]) and set(mask)<=set('01') for mask in masks)
        records.append(line.strip()+' '+' '.join(masks))
    print('N,family,seed,positive_visible,negative_visible,reference_boolean_mismatches,reference_errors,batch_errors')
    with ProcessPoolExecutor(max_workers=4) as pool:
        for row in pool.map(check,records):
            print(*row,sep=',')
            assert row[-3:]==(0,0,0),row
    assert not paired
