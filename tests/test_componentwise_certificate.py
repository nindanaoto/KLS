"""Exact differential audit of the sufficient componentwise certificate."""
import ctypes as C
from fractions import Fraction as F
import math
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='kls-componentwise-test-')
        cls.library = Path(cls.temp.name)/'certificate.so'
        source = Path(__file__).resolve().parents[1]/'src/kls_componentwise.c'
        subprocess.run([os.environ.get('CC','cc'),'-O2','-shared','-fPIC',str(source),'-o',str(cls.library),'-lm'],check=True)
        cls.lib = C.CDLL(str(cls.library))
        cls.check = cls.lib.kls_componentwise_certify
        cls.check.argtypes = [C.c_size_t,C.POINTER(C.c_int64),C.POINTER(C.c_int64),
            C.POINTER(C.c_double),C.POINTER(C.c_double),C.POINTER(C.c_double),
            C.c_int,C.POINTER(C.c_double),C.POINTER(C.c_longdouble)]
        cls.check.restype = C.c_int

    @classmethod
    def tearDownClass(cls): cls.temp.cleanup()

    def compare(self,n,p,rows,a,b,x,transpose=False):
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
        if got==1:self.assertTrue(exact,(p,rows,a,b,x,transpose))
        return got,exact

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


if __name__=='__main__':unittest.main()
