"""Exact differential audit of the sufficient componentwise certificate."""
import ctypes as C
from fractions import Fraction as F
import math
import os
from pathlib import Path
import random
import shlex
import subprocess
import tempfile
import unittest


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='kls-componentwise-test-')
        cls.library = Path(cls.temp.name)/'certificate.so'
        source = Path(__file__).with_name('kls_componentwise_probe.c')
        subprocess.run([os.environ.get('CC','cc'),'-O2',*shlex.split(os.environ.get('CFLAGS','')),
                        '-shared','-fPIC',str(source),'-o',str(cls.library),'-lm'],check=True)
        cls.lib = C.CDLL(str(cls.library))
        cls.check = cls.lib.kls_componentwise_certify
        cls.check.argtypes = [C.c_size_t,C.POINTER(C.c_int64),C.POINTER(C.c_int64),
            C.POINTER(C.c_double),C.POINTER(C.c_double),C.POINTER(C.c_double),
            C.c_int,C.POINTER(C.c_double),C.POINTER(C.c_longdouble)]
        cls.check.restype = C.c_int
        cls.fast = cls.lib.kls_test_fast
        cls.fast.argtypes = cls.check.argtypes
        cls.fast.restype = C.c_int
        cls.bound_text = cls.lib.kls_test_bound
        cls.bound_text.argtypes = [C.POINTER(C.c_longdouble),C.c_char_p,C.c_size_t]
        cls.bound_text.restype = C.c_int
        cls.parallel=cls.lib.kls_test_parallel
        cls.parallel.argtypes=[*cls.check.argtypes,C.c_int]
        cls.parallel.restype=C.c_int

    @classmethod
    def tearDownClass(cls): cls.temp.cleanup()

    def compare(self,n,p,rows,a,b,x,transpose=False,parallel_threads=None):
        pp=(C.c_int64*len(p))(*p);rr=(C.c_int64*len(rows))(*rows)
        aa=(C.c_double*len(a))(*a);bb=(C.c_double*n)(*b);xx=(C.c_double*n)(*x)
        bound=C.c_longdouble()
        got=self.check(n,pp,rr,aa,bb,xx,int(transpose),None,C.byref(bound))
        coeff={}
        for j in range(n):
            for k in range(p[j],p[j+1]):
                i=rows[k];key=(j,i) if transpose else (i,j)
                coeff[key]=coeff.get(key,F(0))+F(a[k])
        residual=list(map(F,b));den=list(map(lambda v:abs(F(v)),b))
        for (i,j),av in coeff.items():
            product=av*F(x[j]);residual[i]-=product;den[i]+=abs(product)
        exact=all(abs(r)<=d/F(10**8) for r,d in zip(residual,den))
        def verify(result,upper):
            if result!=1:return
            self.assertTrue(exact,(p,rows,a,b,x,transpose))
            text=C.create_string_buffer(128)
            length=self.bound_text(C.byref(upper),text,len(text))
            self.assertTrue(0<length<len(text))
            significand,exponent=text.value.decode().lower().split('p')
            whole,point,fraction=significand.partition('.')
            integer=int(whole+fraction,16)
            reported=F(integer)*F(2)**(int(exponent)-4*len(fraction))
            actual=max((abs(r)/d if d else F(0) for r,d in zip(residual,den)),default=F(0))
            self.assertGreaterEqual(reported,actual)
        verify(got,bound)
        fast_bound=C.c_longdouble()
        fast=self.fast(n,pp,rr,aa,bb,xx,int(transpose),None,C.byref(fast_bound))
        verify(fast,fast_bound)
        if fast==1:self.assertEqual(got,1)
        if parallel_threads is not None:
            parallel_bound=C.c_longdouble()
            parallel=self.parallel(n,pp,rr,aa,bb,xx,int(transpose),None,C.byref(parallel_bound),parallel_threads)
            verify(parallel,parallel_bound)
            self.assertEqual(parallel,fast)
            if fast==1:
                left=C.create_string_buffer(128);right=C.create_string_buffer(128)
                self.bound_text(C.byref(fast_bound),left,len(left))
                self.bound_text(C.byref(parallel_bound),right,len(right))
                self.assertEqual(left.value,right.value)
        return got,exact

    def test_parallel_bounds(self):
        if not self.lib.kls_test_parallel_available():self.skipTest('OpenMP fast proof unavailable')
        for threads in (1,2,4,8):
            for n in (1,3,31,65,257):
                values=[math.ldexp(1.+(j%7)/8.,(j%31)-15) for j in range(n)]
                total=sum(map(F,values),F(0));tol=F(1,10**8)
                edge=float(total*(1+tol)/(1-tol))
                for b0 in (float(total),math.nextafter(edge,-math.inf),edge,math.nextafter(edge,math.inf)):
                    for trans in (False,True):
                        p=[0,n]+[n]*(n-1) if trans else list(range(n+1))
                        rows=list(range(n)) if trans else [0]*n
                        self.compare(n,p,rows,values,[b0]+[0.]*(n-1),[1.]*n,trans,threads)

    def test_fast_dispatch(self):
        if not self.lib.kls_test_fast_available():self.skipTest('extended fast proof unavailable')
        pp=(C.c_int64*2)(0,1);rr=(C.c_int64*1)(0)
        aa=(C.c_double*1)(1);bb=(C.c_double*1)(1);xx=(C.c_double*1)(1)
        bound=C.c_longdouble()
        self.assertEqual(self.fast(1,pp,rr,aa,bb,xx,0,None,C.byref(bound)),1)
        # Duplicates must miss the fast proof even when their exact sum is safe.
        dp=(C.c_int64*2)(0,2);dr=(C.c_int64*2)(0,0);da=(C.c_double*2)(.5,.5)
        self.assertEqual(self.fast(1,dp,dr,da,bb,xx,0,None,C.byref(bound)),0)
        self.assertEqual(self.check(1,dp,dr,da,bb,xx,0,None,C.byref(bound)),1)

    def test_random_and_boundaries(self):
        rng=random.Random(5021)
        for _ in range(1500):
            n=rng.randrange(1,7);p=[0];rows=[];a=[]
            x=[math.ldexp(rng.uniform(-1,1),rng.randrange(-500,501)) for i in range(n)]
            b=[math.ldexp(rng.uniform(-1,1),rng.randrange(-500,501)) for i in range(n)]
            for j in range(n):
                for k in range(rng.randrange(1,2*n+1)):
                    rows.append(rng.randrange(n));a.append(math.ldexp(rng.uniform(-1,1),rng.randrange(-500,501)))
                p.append(len(rows))
            self.compare(n,p,rows,a,b,x, bool(rng.randrange(2)))
        for x in [0.,1.,math.nextafter(1.,0.),math.nextafter(1.,2.),
                  1.+2e-8,math.nextafter(1.+2e-8,0.),math.nextafter(1.+2e-8,2.)]:
            self.compare(1,[0,1],[0],[1.],[1.],[x])

    def test_extremes_and_cancellation(self):
        for a in [0.,math.ulp(0.),2.**-1022,1.,2.**1023]:
            got,exact=self.compare(1,[0,1],[0],[a],[a],[1.])
            self.assertTrue(exact);self.assertEqual(got,1)
        got,exact=self.compare(2,[0,1,3],[0,0,1],[1.,2.**100,1.],
                               [1.,1.],[-2.**100,1.])
        self.assertTrue(exact);self.assertEqual(got,1)
        # Inflating |A| by taking abs before duplicate cancellation is unsafe.
        got,exact=self.compare(1,[0,3],[0,0,0],[2.**54,1.,-2.**54],[1.],[2.])
        self.assertFalse(exact);self.assertEqual(got,0)

    def test_near_solutions(self):
        # Exercise the accepting side too, including empty rows, signed
        # duplicates and transpose. Form the RHS exactly before rounding it.
        rng=random.Random(91237)
        accepted=0
        for trans in (False,True):
            for _ in range(300):
                n=rng.randrange(1,9);p=[0];rows=[];a=[]
                x=[math.ldexp(rng.uniform(-1,1),rng.randrange(-20,21)) for i in range(n)]
                exact_b=[F(0) for i in range(n)]
                for j in range(n):
                    for k in range(rng.randrange(0,2*n+1)):
                        i=rng.randrange(n);v=math.ldexp(float(rng.randrange(-16,17)),rng.randrange(-20,21))
                        rows.append(i);a.append(v)
                        dest,src=(j,i) if trans else (i,j)
                        exact_b[dest]+=F(v)*F(x[src])
                    p.append(len(rows))
                b=list(map(float,exact_b))
                got,exact=self.compare(n,p,rows,a,b,x,trans)
                self.assertTrue(exact)
                accepted+=got==1
        # An interval certificate may conservatively reject exact cancellation
        # with a zero denominator; rejection is not evidence of inaccuracy.
        self.assertGreater(accepted,0)

    def test_unique_irregular_rows(self):
        rng=random.Random(29731)
        for trans in (False,True):
            for _ in range(250):
                n=rng.randrange(1,65);p=[0];rows=[];a=[]
                x=[math.ldexp(rng.uniform(-1,1),rng.randrange(-100,101)) for i in range(n)]
                rhs=[F(0) for i in range(n)]
                for j in range(n):
                    selected=rng.sample(range(n),rng.randrange(min(n,12)+1))
                    for i in selected:
                        value=math.ldexp(rng.uniform(-1,1),rng.randrange(-100,101))
                        rows.append(i);a.append(value)
                        dest,src=(j,i) if trans else (i,j)
                        rhs[dest]+=F(value)*F(x[src])
                    p.append(len(rows))
                b=list(map(float,rhs))
                if rng.randrange(2) and b:
                    i=rng.randrange(n);b[i]*=1.+rng.choice((-2e-8,2e-8))
                self.compare(n,p,rows,a,b,x,trans)

    def test_accumulation_boundaries(self):
        # No duplicate coordinates: the fast gate must bound accumulation,
        # not just each individual product, on either side of the tolerance.
        tol=F(1,10**8)
        for n in (1,3,31,128,257):
            values=[math.ldexp(1.+(j%7)/8.,(j%31)-15) for j in range(n)]
            total=sum(map(F,values),F(0))
            edge=float(total*(1+tol)/(1-tol))
            for b0 in (math.nextafter(edge,-math.inf),edge,math.nextafter(edge,math.inf)):
                for trans in (False,True):
                    p=[0,n]+[n]*(n-1) if trans else list(range(n+1))
                    rows=list(range(n)) if trans else [0]*n
                    self.compare(n,p,rows,values,[b0]+[0.]*(n-1),[1.]*n,trans)


if __name__=='__main__':unittest.main()
