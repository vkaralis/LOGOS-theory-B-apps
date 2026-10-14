import tempfile
import unittest
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from logos_b.models import RemModel, WarModel
from logos_b.synthetic import generate
from logos_b.benchmark import canonical

class NumericalContracts(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        generate(self.root,2)
        self.data=json.loads((self.root/'profiles.json').read_text())

    def test_seed_and_bytes(self):
        other=self.root/'repeat';generate(other,2)
        for name in ['profiles.json','observations.csv','mtx_courses.csv']:
            self.assertEqual((self.root/name).read_bytes(),(other/name).read_bytes())

    def test_effect_site_equal_rates(self):
        m=dict(self.data['Rem_models'][0]);m['ke0_min']=m['k1_min'];model=RemModel(m)
        x=np.ones(20)
        sol=solve_ivp(lambda t,y:[-model.k[0]*y[0]+model.w[0]*model.rate,-model.k[1]*y[1]+model.w[1]*model.rate,model.ke*(y[0]+y[1]-y[2])],(0,model.W),[0.,0.,0.],rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(model.values([model.W],x,'ce')[0],sol.y[2,-1],rtol=1e-8)

    def test_effect_site_tail(self):
        model=RemModel(self.data['Rem_models'][0]);x=np.linspace(.5,1.5,20)
        # Resolve both Cp components through independent component ODEs on input bins.
        state=np.zeros(4)
        for j,(a,b) in enumerate(zip(model.edges[:-1],model.edges[1:])):
            sol=solve_ivp(lambda t,y:[-model.k[0]*y[0]+model.w[0]*model.rate*x[j],-model.k[1]*y[1]+model.w[1]*model.rate*x[j],model.ke*(y[0]+y[1]-y[2]),0],(a,b),state,rtol=1e-10,atol=1e-12)
            state=sol.y[:,-1]
        sol=solve_ivp(lambda t,y:[-model.k[0]*y[0],-model.k[1]*y[1],model.ke*(y[0]+y[1]-y[2]),y[2]],(model.W,model.W+2000),state,rtol=1e-10,atol=1e-12)
        self.assertAlmostEqual(model.tail(x),sol.y[3,-1],places=7)

    def test_war_amount_and_absorption(self):
        model=WarModel(self.data['War_models'][0]);t=np.linspace(0,1000,100001)
        cp=model.cp(t,[0,6],[.4,.6])
        auc=np.trapezoid(cp,t)
        self.assertAlmostEqual(auc,model.D*model.scale/model.k,places=3)
        self.assertEqual(float(model.cp([.2],[0],[1])[0]),0.)

    def test_war_turnover_identity(self):
        model=WarModel(self.data['War_models'][0],Rinit=.9)
        times=[0,12,24];x=[.4,.3,.3]
        sol=solve_ivp(lambda t,y:[model.kout*(model.K/(model.K+model.cp([t],times,x)[0])-y[0]),1-y[0]],(0,120),[.9,0],rtol=1e-10,atol=1e-12,max_step=.2)
        self.assertAlmostEqual(model.objective(times,x),sol.y[1,-1],places=6)

    def test_quadratic_reference(self):
        c=canonical()
        self.assertAlmostEqual(c['J'],28.826261662949886,places=7)
        self.assertAlmostEqual(c['B'],12.613903616328,places=9)
        self.assertAlmostEqual(c['delayed_injury_peak_time'],2.01179739054,places=9)
        self.assertGreater(c['J'],c['Jfloor']);self.assertGreater(c['Jfloor'],c['Jref'])

if __name__=='__main__':unittest.main()
