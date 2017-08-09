from py_ibm import *
from plot_trees import plot_tree_pos
import pdb
from tables import *
import json

class Tree(Agent):
    """ Represent individual trees.

    Args:
        position (tuple): (x,y) coordinates of where tree will be stablished.
        world (world object): world where tree will live.
        dbh (float): Diameter at Breast Height (in cm).
        h0 (float): Height-stem diameter relationship parameter.
        h1 (float): Height-stem diameter relationship parameter.
        cl0 (float): Crown Length-Height relationship parameter.
        cd0 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd1 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd2 (float): Crown Diameter-Stem Diameter relationship parameter.
        rho (float): wood density in t of Organic Dry Matter/cubic meter.
        sigma (float): ratio of total aboveground biomass to stem biomass.
        f0 (float): form factor-stem diameter relationship parameter.
        f1 (float): form factor-stem diameter relationship parameter.
        l0 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        l1 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        m (float): Light trasnmission coefficient.
        alpha (float): quantum eciency; initial slope of the type specic light response curve.
        pmax (float): maximum leaf gross photosynthetic rate.
        rg (float): fraction of gross primary production available for growth that is attributed to growth respiration.
        deltaDmax (float): maximum diameter increment.
        DdeltaDmax (float): % of Dmax that reaches deltaDmax.
        age (int): age (years).
        Ftype (str): functional type.
        Mb (float): backgroung mortality probability.
        Dfall (float): mimimum dbh for a tree to fall.
        pfall (float): probabbility that a tree with dbh>Dfall will fall.
        m_max (float): maximum size-dependent mortality for small trees.
        Dmort (float): DBH up to which mortality is increased (for small trees).

    Attributes:

        ID (int):
        Instances (dict):
        DeadTrees (dict):

        patch (tuple): (x,y) patch in which tree is located.
        f (float):
        k (float):
        lday (float):
        phi_act (float):

        H (float): Height (m).
        CL (float): Crown length (m).
        CD (float): Crown diameter (m)
        CA (float): Crown area (sq. m)
        AGB (float): Above ground biomass (tODM)
        LAI (float): Leaf area index.
        lmax (float):
        lmin (float):
        L_mean (float):
        Li (float):
        Mc (float):
        basic_m (float):
        Mb (float):
        Md (float):
        GPP (float): Gross primary production in this time step.
        Rm (float): Respiration maintenance in this timestep.
        position (tuple): (x,y) coordinates of where tree is located.
        world (world object): world where tree lives.
        dbh (float): Diameter at Breast Height (in cm).
        h0 (float): Height-stem diameter relationship parameter.
        h1 (float): Height-stem diameter relationship parameter.
        cl0 (float): Crown Length-Height relationship parameter.
        cd0 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd1 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd2 (float): Crown Diameter-Stem Diameter relationship parameter.
        rho (float): wood density in t of Organic Dry Matter/cubic meter.
        sigma (float): ratio of total aboveground biomass to stem biomass.
        f0 (float): form factor-stem diameter relationship parameter.
        f1 (float): form factor-stem diameter relationship parameter.
        l0 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        l1 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        m (float): Light trasnmission coefficient.
        alpha (float): quantum eciency; initial slope of the type-specic light response curve.
        pmax (float): maximum leaf gross photosynthetic rate.
        rg (float): fraction of gross primary production available for growth that is attributed to growth respiration.
        deltaDmax (float): maximum diameter increment.
        DdeltaDmax (float): % of Dmax that reaches deltaDmax.
        age (int): age (years).
        Ftype (str): functional type.
        Mb (float): backgroung mortality probability.
        max_height (float): maximum height this tree can reach.
        Dfall (float): mimimum dbh for a tree to fall.
        pfall (float): probabbility that a tree with dbh>Dfall will fall.
        m_max (float): maximum size-dependent mortality for small trees.
        Dmort (float): DBH up to which mortality is increased (for small trees).

    """


    ID=0
    Instances={}
    DeadTrees={"FT1":0,"FT2":0,"FT3":0,"FT4":0,"FT5":0,"FT6":0}

    @classmethod
    def TreesAboveDBH(self,dbh):
        """Calculate the total number of trees with DBH above the speciefied threshold.

        Returns:
            (int): Number of trees.
        """

        trees=[k for k in self.Instances.keys() if self.Instances.get(k).DBH>=dbh]
        return len(trees)

    @classmethod
    def Total_AGB(self):
        """ Calculate the total Above Ground Biomass for all living trees.

        Returns:
            (float): total AGB in t of Organic Dry matter (tODM).
        """

        Tagb=0
        trees=[k for k in self.Instances.keys()]
        for i in trees:
            t=Tree.Instances.get(i)
            Tagb+=t.AGB
        return Tagb

    @classmethod
    def Cgpp(self):
        """ Calculate the total Carbon absorbed as Gross Primary Production.

        Returns:
            (float): Total Carbon absorbed.
        """


        Cgpp=0
        trees=[k for k in self.Instances.keys()]
        for i in trees:
            t=Tree.Instances.get(i)
            Cgpp+=t.GPP
        return 0.44*Cgpp

    @classmethod
    def Cr(self):
        """ Calculate the total Carbon released through respiration.

        Returns:
            (float): Total Carbon released by all living trees.
        """
        Cr=0
        trees=[k for k in self.Instances.keys()]
        for i in trees:
            t=Tree.Instances.get(i)
            Cr+=(t.Rm+t.rg*(t.GPP-t.Rm))
        return 0.44*Cr

    @classmethod
    def Mortality(self):
        """
        Kill trees by calling their stochastic_M() method.
        """
        trees=[k for k in Tree.Instances.keys()]
        for t in trees:
             if Tree.Instances.get(t).stochastic_M(): Tree.Instances.get(t).remove_me()


    @classmethod
    def BackgroundMortality(self):
        """
        Kill trees by calling their stochastic_Bm() method.
        """
        trees=[k for k in self.Instances.keys()]
        for t in trees:
             if Tree.Instances.get(t).stochastic_Bm(): Tree.Instances.get(t).remove_me()

    @classmethod
    def CrowdingMortality(self):
         """
         Kill trees by calling their stochastic_Cm() method.
         """
         trees=[k for k in self.Instances.keys()]
         for t in trees:
              if Tree.Instances.get(t).stochastic_Cm(): Tree.Instances.get(t).remove_me()

    @classmethod
    def DamageMortality(self):
          """
          Kill trees by calling their stochastic_Dm() method.
          """
          trees=[k for k in self.Instances.keys()]
          for t in trees:
               if Tree.Instances.get(t).stochastic_Dm(): Tree.Instances.get(t).remove_me()


    @classmethod
    def ActiveAgents(self):
        return [ag.id for ag in self.Instances.values() if ag.active == True]

    @classmethod
    def TreesPerPatch(self):
        """Identify which trees are in each patch.

        Returns:
            (dict): keys are patch coordinates '(x,y)' and values are lists of trees ids in that patch.

        """
        patches={}
        for t in self.Instances.values():
            if patches.get(t.patch)==None:
                patches[t.patch]=[t.id]
            else:
                patches[t.patch].append(t.id)
        return patches

    @classmethod
    def MaxHeight(self):
        """ Finds the height of the tallest tree in the current time step.

        Returns:
            (float): height of the tallest tree.

        """
        return max([t.H for t in self.Instances.values()],default=0)


    def __repr__(self):
        return "Tree type: {0}".format(self.Ftype)

    def __init__(self,position,world,dbh,wood_density,h0,h1,cl0,cd0,cd1,cd2,rho,sigma,f0,f1,l0,l1,m,alpha,pmax,max_agb,rg,deltaDmax,DdeltaDmax,age,Ftype,Mb,Dfall=45, pfall=0.3,m_max=0.12,Dmort=10,max_dbh=None,id=None):

        super().__init__(position,world,id=id)
        #self.AddInstance(self.id)
        #Tree.Instances[Ftype+' '+str(self.id)]=self

        self.patch=(int(np.floor(self.position[0])),int(np.floor(self.position[1])))
        self.DBH=dbh
        self.age=age
        self.Ftype=Ftype
        self.wood_density=wood_density
        #type specific parameter used to calculate tree height
        self.h0=h0
        self.h1=h1
        #type specific parameter used to calculate crown length
        self.cl0=cl0
        #type specific parameter used to calculate crown diameter
        self.cd0=cd0
        self.cd1=cd1
        self.cd2=cd2
        #type specific parameter used to calculate aboveground biomass
        self.rho=rho
        self.sigma=sigma
        self.f0=f0
        self.f1=f1
        self.f=self.calculate_f()
        #type specific parameter used to calculate the leaf area index
        self.l0=l0
        self.l1=l1
        self.k=self.world.topology.k
        self.lday=self.world.topology.lday
        self.phi_act=self.world.topology.phi_act
        self.m=m
        self.alpha=alpha
        self.pmax=pmax
        self.m_max=m_max
        self.Dmort=Dmort
        self.Dfall=Dfall
        self.pfall=pfall
        self.Dmax=max_dbh/100.
        self.deltaDmax=deltaDmax
        self.DdeltaDmax=DdeltaDmax
        self.alpha0=self.calculate_alpha0()
        self.alpha1=self.calculate_alpha1()
        self.rg=rg

        self.H=0
        self.CL=0
        self.CD=0
        self.CA=0
        self.AGB=0
        self.LAI=0
        self.lmax=0
        self.lmin=0
        self.L_mean=0
        self.Li=0
        self.Mc=0
        self.basic_m=Mb
        self.Mb=0
        self.Md=0
        self.GPP=0 #Gpp in this time step
        self.Rm=0 #Rm in this timestep
        #self.update()
        self.world.topology.trees_per_patch[self.patch].append(self.id)


    def update(self):
        """ Updates geometry related parameters and mortality probabilities.
        """
        self.H=self.H_from_DBH()
        self.CL=self.CL_from_H()
        self.CD=self.CD_from_DBH()
        self.CA=self.CA_from_CD()
        if self.age==0:
            self.AGB=self.AGB_from_DBH(D=self.DBH/100)
        self.LAI=self.LAI_from_DBH()

        self.lmax=self.H/self.world.topology.delta_h #eq.31
        self.lmin=(self.H-self.CL)/self.world.topology.delta_h #eq.32
        self.L_mean=(self.LAI*self.CA)/(self.lmax-self.lmin if self.lmax-self.lmin!=0 else 0.5 )#eq.34
        self.Li=self.Lind()
        self.Mc=self.calculate_Rc()
        self.Mb=max(0,self.basic_m+self.calculate_Ms())
        self.Md=0


        attributes={"DBH":self.DBH,
                "H":self.H,
                "CL":self.CL,
                "CD":self.CD,
                "CA":self.CA,
                "AGB":self.AGB,
                "LAI":self.LAI,
                "lmax":self.lmax,
                "lmin":self.lmin,
                "L_mean":self.L_mean,
                "Li":self.Li,
                "Mc":self.Mc,
                "Mb":self.Mb,
                "age":self.age,}

        return (attributes)

    def import_attributes(self,attributes):
        self.DBH=attributes["DBH"]
        self.H=attributes["H"]
        self.CL=attributes["CL"]
        self.CD=attributes["CD"]
        self.CA=attributes["CA"]
        self.AGB=attributes["AGB"]
        self.LAI=attributes["LAI"]
        self.lmax=attributes["lmax"]
        self.lmin=attributes["lmin"]
        self.L_mean=attributes["L_mean"]
        self.Li=attributes["Li"]
        self.Mc=attributes["Mc"]
        self.Mb=attributes["Mb"]
        self.age=attributes["age"]
        self.Md=0
        #self.age+=1





    def H_from_DBH(self):
        """Calculate tree height calculated according to eq.1 in SI-Fisher et al.(2015).

        Returns:
            (float): H.
        """
        return self.h0*self.DBH**self.h1

    def CL_from_H(self):
        """Calculate tree crown_length (CL) calculated according to eq.2 in SI-Fisher et al.(2015)

        Returns:
            (float): CL.
        """
        return self.cl0*self.H

    def CD_from_DBH(self):
        """Calculate tree crown diameter (CD) calculated according to eq.3 in SI-Fisher et al.(2015).

        Returns:
            (float): CD.
        """
        return (self.cd0*self.DBH**self.cd1)-self.cd2

    def CA_from_CD(self):
        """Calculate tree crown area (CA) calculated according to eq.4 in SI-Fisher et al.(2015)

        Returns:
            (float): CA.
        """
        return (np.pi/4)*(self.CD**2)

    def calculate_volume(self):
        """ Calculate the volume of the stem based on DBH.

        Returns:
            (float): volume.
        """

        return self.H*(np.pi*((self.DBH/100)/2)**2)


    def AGB_from_DBH(self,D):
        """ Calculate tree aboveground biomass (AGB) calculated according to eq.5 in SI-Fisher et al.(2015)

        Returns:
            (float): AGB.
        """

        return (np.pi/4)*(D**2)*self.H*self.f*(self.rho/self.sigma)

    def calculate_f(self):
        """Calculate the form factor (f) calculated according to eq.6 in SI-Fisher et al.(2015)

        Returns:
            (float): f.

        """

        return self.f0*self.DBH**self.f1

    def LAI_from_DBH(self):
        """Calculate the leaf area index (LAI) calculated according to eq.7 in SI-Fisher et al.(2015)

        Returns:
            (float): LAI.
        """

        return self.l0*self.DBH**self.l1

    def Lind(self):
        """Calculate the incoming radiation on top of lmax layer the tree is reaching according to eq.36 in SI-Fisher et al.(2015).

        Returns:
            (float): Incoming radiation on top of this tree.
        """

        #crown_layers=(np.ceil(self.lmax)-np.floor(self.lmin))/self.world.topology.delta_h
        patch_LI=self.world.topology.LAIs[self.patch]
        top_layer=int(np.ceil(self.lmax))#/self.world.topology.delta_h)
        patch_LI=patch_LI[top_layer:]

        #self.patch_LI=sum(patch_LI)

        return self.world.topology.I0*np.exp(-self.k*sum(patch_LI))


#    def Lleaf(self,L):
#        """
#        Return the incoming radiation on top of leaves in layerL, according to eq.38 in SI-Fisher et al.(2015)
#        """

#        return (self.k/(1-self.m))*self.Li*np.exp(-self.k*L)

#    def Pleaf(self,L):
#        """
#        Return the gross photosynthetic rate on top of leaves in layer L, according to eq.37 in SI-Fisher et al.(2015)
#        """
#        lleaf=self.Lleaf(L)

#        return (self.alpha*lleaf*self.pmax)/(self.alpha*lleaf+pmax)

    def Pind(self,Li):
        """ Calculate the interim gross photosynthetic rate of one tree per year, according to eq.40 in SI-Fisher et al.(2015)

        Returns:
            (float) Gross photosynthetic rate.
        """

        pmk=self.pmax/self.k
        aki=self.alpha*self.k*Li
        pm=self.pmax*(1-self.m)

        pind=pmk*np.log((aki+pm)/(aki*np.exp(-self.k*self.LAI)+pm))

        #Tons of organic dry meter per year (todm.y-1)
        return pind*self.CA*60*60*self.lday*self.phi_act*np.power(10.0,-12)*2.27*44


    def calculate_rm(self):
            """Calculate the maintenance respiration rate(rm) according to eq.46 in SI-Fisher et al.(2015).

            Returns:
                (float): rm.
            """

            max_Li=self.world.topology.I0

            rm=(1/self.AGB)*(self.Pind(max_Li)-(max(0,(self.AGB_from_DBH(D=(self.DBH/100)+self.growth(D=self.DBH/100.))-self.AGB))/(1-self.rg)))

            return rm

    def growth(self,D):
        """Calculate the yearly diameter growth according to eq.48 in SI-Fisher et al.(2015). Assumes full availabilty of resources

        Returns:
            (float):diameter growth.
        """

        return (self.alpha0*D*(1-(D/self.Dmax)))*np.exp(-self.alpha1*D)

    def calculate_alpha0(self):
        """Calculate the alpha0 growth paramenter, according equations described in sctions F-5 (page 24) of SI-Fisher et al.(2015).

        Returns:
            (float): alpha 0.

        """


        exp=np.exp((self.Dmax-2*(self.DdeltaDmax*self.Dmax))/(self.Dmax-(self.DdeltaDmax*self.Dmax)))
        alpha0=exp*self.Dmax*self.deltaDmax
        alpha0=alpha0/((self.Dmax-(self.DdeltaDmax*self.Dmax))*(self.DdeltaDmax*self.Dmax))

        return alpha0

    def calculate_alpha1(self):
        """Calculate alpha1 growth paramenter, according equations described in sctions F-5 (page 24) of SI-Fisher et al.(2015).

        Returns:
            (float): alpha1.
        """

        alpha1=self.Dmax-2*(self.DdeltaDmax*self.Dmax)
        alpha1=alpha1/(self.Dmax*(self.DdeltaDmax*self.Dmax)-(self.DdeltaDmax*self.Dmax)**2)

        return alpha1

    def Biomass_gain(self):
        """ Calculate Aboveground Biomass (AGB) gain according to eq.43 in SI-Fisher et al.(2015).

        Returns:
            (float): AGB.
        """

        self.GPP=self.Pind(self.Li)
        self.Rm=self.calculate_rm()*self.AGB
        #self.Rm=self.r0*self.AGB+self.r1*self.AGB**2+self.r2*self.AGB**3
        gain=(1-self.rg)*(self.GPP-self.Rm)
        self.AGB+=gain
        return gain

    def DBH_from_Biomass(self,B):
        """ Calculate diameter at breast height from above ground biomass.

        Returns:
            (float): the DBH calculated from AGB.
        """


        p=(np.pi/4)*self.H*self.f*(self.rho/self.sigma)
        #if p>0 and B>0:
        dbh=np.sqrt(B/p)
        if isnan(dbh) or dbh<=0: return 0 #if Biomass is negative, new DBH will be 0
        return dbh

    def increase_DBH(self):
        """Increase tree DBH according to gain in biomass.

        Note: If DBH decreases to zero, remove tree.

        Returns:
            None.
        """
        gain=self.Biomass_gain()
        if round(self.AGB,3) <=0 or gain<0:
            self.remove_me()
        else:
            #new_B=self.AGB+self.Biomass_gain()
            self.DBH=self.DBH_from_Biomass(self.AGB)*100
            if round(self.DBH,3)<=2: self.remove_me()
        #if self.age>150: self.remove_me()

        #new_B=self.AGB+self.Biomass_gain()
        #self.DBH=self.DBH_from_Biomass(new_B)*100
        #if self.DBH<=0: self.remove_me()

    def calculate_Rc(self):
        patch_CCA=self.world.topology.CCAs[self.patch]
        max_CCA=max(patch_CCA[int(np.ceil(self.lmin)):int(np.ceil(self.lmax))],default=self.CA)

        if max_CCA==0: return 1.0
        rc=1/max_CCA
        if rc>=1.0:
            return 0
        else:
            return (1-rc)

    def calculate_Ms(self):
        """Calculate the size-dependent mortality for small trees according to model description in Ruger et al. 2007. (SI, page 4).

        Returns:
            (float): Ms.
        """
        Ms=0


        if self.DBH <self.Dmort:
            Ms=self.m_max-self.m_max*(self.DBH/self.Dmort)
        #elif self.DBH>=self.max_dbh:
        #    Ms=self.m_max

        return Ms


    def stochastic_M(self):
        """Stochastic mortality.

        Returns:
            (bools): True if tree should be removed.
        """
        #if self.DBH>=self.max_dbh:
        #    return True #np.random.random()<= self.m_max
        w_Mb=0
        w_Mc=1
        w_Md=1


        M=((w_Mc*self.Mc+w_Mb*self.Mb+w_Md*self.Md)/(w_Mc+w_Mb+w_Md))#
        return np.random.random()<=M#*self.world.mortality_factor

    def stochastic_Cm(self):
        """Stochastic crownding mortality.

        Returns:
            (bool): True if tree dies from crowding.
        """
        return np.random.random()<=min(self.Mc,1.0)

    def stochastic_Bm(self):
        """Stochastic Background mortality.

        Returns:
            (bool): True if tree dies from background mortality.
        """
        return np.random.random()<=min(self.Mb,1.0)

    def stochastic_Dm(self):
        """Stochastic damage mortality.

        Returns:
            (bool): True if tree dies from damage.
        """
        return np.random.random()<=min(self.Md,1.0)

    def log_me(self):
        """Cut down this tree.

            Randomly determines where tree is going to fall and inflict damage in the the trees hit by the crown. Only trees with DBH< 50 cm are damaged.

            The carbon of this trees is not added to any Stock.
        """


        subclass="Tree_"+self.Ftype
        globals()[subclass].Indices.remove(self.id)
        self.world.topology.trees_per_patch[self.patch].remove(self.id)
        #self.world.Smort+=self.AGB*0.44

        topology=self.world.topology

        dir=np.random.randint(1,360)
        x=self.position[0]+self.H*np.sin(np.deg2rad(2*np.pi*(dir/360)))
        y=self.position[1]+self.H*np.cos(np.deg2rad(2*np.pi*(dir/360)))
        md=min(self.CA/topology.patch_area,1)

        target_patch=np.floor((x,y))
        target_patch=(int(target_patch[0]),int(target_patch[1]))
        target_patch=topology.in_bounds(target_patch)



        neighborhood=topology.hood(cell=target_patch,remove_center=False)
        trees_in_neighborhood=topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if topology.point_within_circle(x,y,(self.CD/2),p[0],p[1],topology.scaled_distance,scale=topology.patch_area**0.5)]

        for t in trees_in_circle:
            if Tree.Instances.get(t).DBH <=50:
                 Tree.Instances.get(t).Md=md



        del Tree.Instances[self.id]


    def remove_me(self):
        """
        Kill this tree.
        """
        subclass="Tree_"+self.Ftype
        globals()[subclass].Indices.remove(self.id)
        if self.DBH>self.Dfall and np.random.random()<=self.pfall:
            self.fall()
        self.world.topology.trees_per_patch[self.patch].remove(self.id)
        self.world.Smort+=max(0,self.AGB)*0.44
        Tree.DeadTrees[self.Ftype]+=1
        del Tree.Instances[self.id]

    def fall(self):
        """
        Inflict damage on other trees affected by the falling one.
        """
        topology=self.world.topology

        dir=np.random.randint(1,360)
        x=self.position[0]+self.H*np.sin(np.deg2rad(2*np.pi*(dir/360)))
        y=self.position[1]+self.H*np.cos(np.deg2rad(2*np.pi*(dir/360)))
        md=min(self.CA/topology.patch_area,1)

        #Damage caused on trees near the falling trees, due to lianas

        neighborhood=topology.hood(cell=self.patch,remove_center=False)
        trees_in_neighborhood=topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if topology.point_within_circle(self.position[0],self.position[1],(self.CD/2),p[0],p[1],topology.scaled_distance,scale=topology.patch_area**0.5)]

        for t in trees_in_circle:
            if Tree.Instances.get(t).H <=(self.H+1):
                 Tree.Instances.get(t).Md=md*(topology.surface[topology.dim_names["LianaCoverage"]][self.patch])



        target_patch=np.floor((x,y))
        target_patch=(int(target_patch[0]),int(target_patch[1]))
        target_patch=topology.in_bounds(target_patch)



        neighborhood=topology.hood(cell=target_patch,remove_center=False)
        trees_in_neighborhood=topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if topology.point_within_circle(x,y,(self.CD/2),p[0],p[1],topology.scaled_distance,scale=topology.patch_area**0.5)]

        for t in trees_in_circle:
            if Tree.Instances.get(t).H <=(self.H+1):
                 Tree.Instances.get(t).Md=md
