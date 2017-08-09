from py_ibm import *
from plot_trees import plot_tree_pos
import pdb
from tables import *
import json

class Tree_World(World):
    """ Provide observer methods and control the schedule and execution of simulations.


    Args:
        topology (Tree_Grid instance): grid object on which agents will live.
        AET (int): actual evapotranspiration in the previous year in mm.
        Assumed to be constant every year.

    Attributes:
        Smort (float): Carbon of all trees that died in one year [tC]
        Sdead (float): Dead wood carbon stock [tC]
        tsdead (float): Anual decomposition rate
        tsdead_A (float): Proportion of decomposed dead wood that is released to the atmosphere.
        Sslow (float): Slow decomposition portion of soil stock [tC]
        tsdead_Sslow (float): Proportion of decomposed dead wood that is tranferred to the slow decomposition Soil stock
        self.tsslow_A (float):Proportion of slow decomposition soil stock that is released to the atmosphere.
        Sfast (float):Fast decomposition portion of soil stock [tC]
        tsfast_A (float):Proportion of fast decomposition soil stock that is released to the atmosphere.
        tsdead_Sfast (float): Proportion of decomposed dead wood that is tranferred to the fast decomposition Soil stock.
        Cgpp (float): Carbon captured in the gross primary productivity of the living forest
        Cr (float):Carbon released by the total respiration of the living forest (i.e. for maintenance and growth)

    Usage:
        >>>grid01=Tree_Grid(x_max=5,y_max=5,delta_h=0.5,patch_area=400,I0=860,lday=12,phi_act=360,k=0.7)
        >>>world01=Tree_World(topology=grid01)
        >>> print(world01)
        World object
    """

    def __init__(self, topology, AET=1300):
        super().__init__(topology)
        self.mortality_factor=1
        self.AET=AET
        self.Smort=0
        self.Sdead={"current":0,"previous":0}
        self.tsdead=min(1.0,np.power(10,(-1.4553+0.0014175*self.AET))/12)
        self.recruit_factor=1.0
        self.dwood_emission_p=0.7
        self.tsdead_A=self.dwood_emission_p*self.tsdead
        self.Sslow={"current":0,"previous":0}
        self.tsdead_Sslow=0.015*0.3*self.tsdead
        self.tsslow_A=1/750
        self.Sfast={"current":0,"previous":0}
        self.tsfast_A=1/15
        self.tsdead_Sfast=0.985*0.3*self.tsdead
        self.smort_log=[]
        self.Cgpp=0
        self.Cr=0
        self.NEE=[]
        self.Stocks={"Sslow":[0,],
                    "Sfast":[0,],
                    "Dwood":[0,],
                    "AGB":[0,],}
        self.tree_subclasses={"FT1":Tree_FT1,
                        "FT2":Tree_FT2,
                        "FT3":Tree_FT3,
                        "FT4":Tree_FT4,
                        "FT5":Tree_FT5,
                        "FT6":Tree_FT6}



    def run_simulation(self,n,logging_years,min_dbh,max_dbh,vol,h5file=None,plot_func=None,**kwargs):
        """Run the code in run_schedule() 'n' times.

        Args:
            n (int): Number of simulations to be run.
            logging_years (range): The years in which selective logging occurs.
            min_dbh (float): minimum dbh for a tree to be considered suitable for logging.
            max_dbh (float): maximum dbh for a tree to be considered suitable for logging.
            vol (float): total volume of timber to be extracted by logging event.
            h5file (str): path to the HDF5 file in which model outputs will be saved
            plot_func (function): a function that will be run every step. Usually a plotting function that saves a image file to the disk. Must accept a 'step' argument,which will receive the step number and can be used to generate the file name.
            **kwargs: List of keyword arguments to be passed to plot_func in addition to 'step'.

        Returns:
            None
        """
        if h5file:

            database=open_file(h5file,"r+")

            Pop_table=database.root.sim_1.sys_lvl.Pop
            Stocks_table=database.root.sim_1.sys_lvl.Stocks
            NEE_table=database.root.sim_1.sys_lvl.NEE
            #Emissions_table=database.root.sim_1.sys_lvl.Emissions


            pop_r=Pop_table.row
            stocks_r=Stocks_table.row
            nee_r=NEE_table.row
            #emissions_r=Emissions_table.row

        self.topology.update()
        for t in Tree.Instances.values():
            t.update()
        self.stablish_new_trees()

        for i in range(n):

            print("step=",i,"N=",len(Tree.Instances),"DBH>10",Tree.TreesAboveDBH(10.0), "DTrees", sum(Tree.DeadTrees.values()))
            if h5file:
               pop_r['FT1']=len(Tree_FT1.Indices)
               pop_r['FT2']=len(Tree_FT2.Indices)
               pop_r['FT3']=len(Tree_FT3.Indices)
               pop_r['FT4']=len(Tree_FT4.Indices)
               pop_r['FT5']=len(Tree_FT5.Indices)
               pop_r['FT6']=len(Tree_FT6.Indices)
               pop_r.append()

               stocks_r['agb']=Tree.Total_AGB()*.44
               stocks_r['dead_wood']=self.Sdead["current"]
               stocks_r['soil_slow']=self.Sslow["current"]
               stocks_r['soil_fast']=self.Sfast["current"]
               stocks_r.append()



            self.topology.update()
            for t in Tree.Instances.values():
                t.update()
            self.topology.update()


            if h5file:

               nee_r['dead_wood']=self.tsdead_A*self.Sdead["current"]
               nee_r['soil_slow']=self.tsslow_A*self.Sslow["current"]
               nee_r['soil_fast']=self.tsfast_A*self.Sfast["current"]
               nee_r['gpp']=self.Cgpp
               nee_r['living_trees']=self.Cr
               nee_r['smort']=self.Smort
               nee_r['t_deadwood_sslow']=self.tsdead_Sslow*self.Sdead["previous"]
               nee_r['t_sslow_sfast']=self.tsdead_Sfast*self.Sdead["previous"]



               nee=self.run_schedule(step=i,n=n,h5_database=database)

               nee_r['total_emissions']=-nee+self.Cgpp

               nee_r['nee']=nee
               nee_r.append()
               NEE_table.flush()
               Pop_table.flush()
               Stocks_table.flush()


            else:
                self.NEE.append(self.run_schedule(step=i,n=n))#,h5_database=database)
                self.Stocks["Sslow"].append(self.Sslow["current"])
                self.Stocks["Sfast"].append(self.Sfast["current"])
                self.Stocks["Dwood"].append(self.Sdead["current"])
                self.Stocks["AGB"].append(Tree.Total_AGB()*.44)


            if not plot_func==None: plot_func(step=i,**kwargs)

            if i in logging_years:
                ps=list(self.topology.trees_per_patch.keys())
                st=self.suitable_trees(ps,min_dbh,max_dbh)
                self.log_trees(st,vol)


            self.stablish_new_trees()
        if h5file: database.close()

    def run_schedule(self,step,n, h5_database=None):
        """Schedule of actions to be run in every time step

        Args:
            step (int): The current step. Received from run_simulation().
            h5_database (str): path to the HDF5 file in which model outputs will be saved

        Returns:
            (float) the Net Ecosystem Exchange (as calculater by 'calculate_NEE') for thar step.
        """
        #self.topology.update()

        list_of_patches=self.topology.trees_per_patch.keys()
        cohorts=self.topology.cohorts(list_of_patches)
        FTs=list(cohorts.keys())


        for ft in FTs:
            ages=list(cohorts[ft].keys())
            for age in ages:
                base_tree_index=cohorts[ft][age][0]
                base_tree=Tree.Instances.get(base_tree_index)
                base_tree.increase_DBH()
                base_tree.age+=1
                attributes=base_tree.update()

                for tree_index in cohorts[ft][age][1:]:
                    tree=Tree.Instances.get(tree_index)
                    tree.import_attributes(attributes)



        #print("DBH:",Tree.Instances[125].DBH,"AGB:",Tree.Instances[125].AGB )

        #trees=[k for k in Tree.Instances.keys()]
        #for i in trees:
        #    t=Tree.Instances.get(i)
            #print("dbh:",t.id)
        #    t.increase_DBH()

        #

        if h5_database:
            Ind_table=h5_database.root.sim_1.ind_lvl.Ind
            ind_r=Ind_table.row
            trees=[k for k in Tree.Instances.keys()]
            for i in trees:
               t=Tree.Instances.get(i)
            #    #print("up:",t.id)
            #    t.update()
            #    t.age+=1


               ind_r["step"]=step
               ind_r["ind_id"]=t.id
               ind_r["pos_x"]=round(t.position[0],2)
               ind_r["pos_y"]=round(t.position[1],2)
               ind_r["age"]=t.age
               ind_r["dbh"]=t.DBH
               ind_r["agb"]=t.AGB
               ind_r["height"]=t.H
               ind_r["crown_area"]=t.CA
               ind_r["FT"]=t.Ftype
               ind_r.append()

            Ind_table.flush()

        self.Cgpp=Tree.Cgpp()
        self.Cr=Tree.Cr()
        Tree.Mortality()
        #Tree.BackgroundMortality()
        # Tree.CrowdingMortality()
        # Tree.DamageMortality()
        return self.calculate_NEE()


    def create_agents(self,FT,n,pos=None, ids=None, DBHs=None,ages=None, **kwargs):
        """ Creates 'n' trees of type 'FT'. A list of positions, ids and DBHs can be provided (Useful to recreate trees from a previous model run, for example).

        Note: overwrites methods create_agents from class Agent.

        Args:
            agent_type (class): the type (class) of agents to be created
            n (int): number of agents to be created
            pos (list,None): a list of 'n' tuples '(x,y)' with the coordinates in which agents will be created. If set to None, agents are created in random positions.
            ids (list): a list of 'n' ints to be used as tree ids.
            DBHs: a list of 'n' floats to be used as tree DBHs.
            kwargs: the arguments to be passed to the 'agent_type' class. All agents will be created with this same set of arguments.

        Returns:
            None.
        """
        if pos==None:
            pos=[(np.random.rand()*self.topology.x_max,
            np.random.rand()*self.topology.y_max)
            for i in range(n)]

        if ids==None:
            for i in range(n):
                FT(position=pos[i],**kwargs)

        else:
            for i in range(n):
                FT(position=pos[i],dbh=DBHs[i], id=ids[i],age=ages[i],**kwargs)


    def seedbank_from_file(self, input_file):
        with open(input_file,'r') as json_file:
            self.topology.seedbank=json.load(json_file)



    def seedbank_to_file(self,output_file):
        with open(output_file,'w') as json_file:
            json.dump(self.topology.seedbank,json_file)


    def model_status_from_file(self, input_file):
        with open(input_file,'r') as json_file:
            trees=json.load(json_file)

        for ft in self.tree_subclasses.keys():
            ages=[]
            positions=[]
            ids=[]
            DBHs=[]
            for tree in trees[ft]:
                positions.append(tuple(tree['position']))
                ids.append(tree['id'])
                DBHs.append(tree['DBH'])
                ages.append(tree['age'])

            self.create_agents(self.tree_subclasses[ft],len(trees[ft]),pos=positions,DBHs=DBHs,ids=ids, ages=ages, world=self)

        self.topology.update()
        for t in Tree.Instances.values():
            t.update()
            t.AGB=t.AGB_from_DBH(t.DBH/100)
        Tree.ID=max(Tree.Instances.keys())+1


    def model_status_to_file(self,output_file):
        trees_per_type={k:[] for k in self.tree_subclasses.keys()}
        for t_id,t in Tree.Instances.items():
            trees_per_type[t.Ftype].append({"id":t_id, "position":t.position,"DBH":t.DBH,"age":t.age})

        with open(output_file,'w') as json_file:
            json.dump(trees_per_type,json_file)


    def suitable_trees(self, patches, min_dbh, max_dbh):
        """Search for suitable trees for logging in the given patches.

        Args:
            patches (list): A list of patches (tuples in the format (x,y))
            min_dbh (float): minimum dbh for a tree to be considered suitable for logging.
            max_dbh (float): max dbh for a tree to be considered suitable for logging.

        Returns:
            A dictionary with patches as keys ((x,y)) and lists of suitable tree ids as values
        """

        suitables_per_patch={}
        for p in patches:
            suitables_per_patch[p]=[t for t  in self.topology.trees_per_patch[p] if Tree.Instances[t].DBH>=min_dbh and Tree.Instances[t].DBH<=max_dbh]

        suitables_per_patch={k:v for k,v in suitables_per_patch.items() if v}
        return suitables_per_patch

    def log_trees(self, suitable_trees, total_vol):
        """Randomly log 'suitable_trees' until the 'total_vol' is reached.

        Args:
            suitable_trees (dict): A dictionary with patches as keys ((x,y)) and lists of suitable tree ids as values (resulting from suitable_trees())
            total_vol (float): the total vol of timber to be logged

        Returns:
            None
        """

        total_logged=0
        while total_logged<total_vol and len(suitable_trees):
            ps=list(suitable_trees.keys())
            p=ps[np.random.randint(len(ps))]
            vols={t:Tree.Instances[t].calculate_volume() for t in suitable_trees[p]}
            #pdb.set_trace()
            biggest_tree=max(vols, key=lambda x: vols[x])
            total_logged+=vols[biggest_tree]
            Tree.Instances[biggest_tree].log_me()
            suitable_trees[p].remove(biggest_tree)
            if not suitable_trees[p]:
                del suitable_trees[p]



    def define_seedbank(self,Nseed):
        """Randomly determines how many seed each patch will receive.

        Args:
            Nseed (int):The total number of seed to be distributed among all patches.

        Returns:
            A dictionary in which keys are patches ('(x,y)') and values are the number of seeds
        """

        seedbank={}
        for s in range(Nseed):
            p=(np.random.randint(self.topology.x_max),np.random.randint(self.topology.y_max))
            if seedbank.get(p)==None:
                seedbank[p]=1
            else:
                seedbank[p]+=1

        return seedbank

    def seeds_pos(self,est_seeds):
        """Set a random position within the patch for each seed in est_seeds.

        Args:
            est_seed (dict):A dictionary with patches as keys ('(x,y)') and the number of established seeds as values.

        Returns:
            A list of seed positions within patches.
        """

        pos=[]
        for i in est_seeds.items():
            for j in range(i[1]):
                seed_pos=i[0]+np.random.rand(2)
                pos.append(self.topology.in_bounds(seed_pos))
        return pos



    def germinate_suitable_seeds(self,FT,pos):
        self.create_agents(self.tree_subclasses[FT],len(pos),pos=pos, world=self, dbh=2)


    def stablish_new_trees(self):
        """Create new tree individuals. Use 'define_seedbank()', 'topology.seed_establishment' and 'seeds_pos'.

            Note: This method can be substitued by an equivalent in order to incorporate more realistic representations of seed dispersal.

        Returns:
            None.
        """


        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT1.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT1.Iseed,h0=Tree_FT1.h0,h1=Tree_FT1.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT1"]=pos
        self.germinate_suitable_seeds(FT="FT1",pos=pos)
        #self.create_agents(Tree_FT1,n,pos=pos, world=self,dbh=2)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT2.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT2.Iseed,h0=Tree_FT2.h0,h1=Tree_FT2.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT2"]=pos
        self.germinate_suitable_seeds(FT="FT2",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT3.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT3.Iseed,h0=Tree_FT3.h0,h1=Tree_FT3.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT3"]=pos
        self.germinate_suitable_seeds(FT="FT3",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT4.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT4.Iseed,h0=Tree_FT4.h0,h1=Tree_FT4.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT4"]=pos
        self.germinate_suitable_seeds(FT="FT4",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT5.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT5.Iseed,h0=Tree_FT5.h0,h1=Tree_FT5.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT5"]=pos
        self.germinate_suitable_seeds(FT="FT5",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT6.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT6.Iseed,h0=Tree_FT6.h0,h1=Tree_FT6.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT6"]=pos
        self.germinate_suitable_seeds(FT="FT6",pos=pos)





    def calculate_Sdead(self):
        """ Calculate the amount of carbon (tC) in the dead wood stock based
            on the trees that died in the previous year (Smort).

            Returns:
                None.

        """
        #self.Smort/=2
        self.Sdead["current"]=max(0,self.Sdead["previous"]+(self.Smort)-self.tsdead*self.Sdead["previous"])
        self.Smort=0

    def calculate_Sslow(self):
        """ Calculate the amount of carbon (tC) in the slow decomposition soil
            stock.

            Calculation is based on decomposition rates for the dead wood (tsdead_Sslow) and the rate with which carbon is transferred from this stock to the atmosphere (tsslow_A).

            Returns:
                None.
        """

        self.Sslow["current"]=self.Sslow["previous"]+self.tsdead_Sslow*self.Sdead["previous"]-self.tsslow_A*self.Sslow["previous"]

    def calculate_Sfast(self):
        """ Calculate the amount of carbon (tC) in the fast decomposition soil
            stock.

            Calculation is based on decomposition rates for the dead wood (tsdead_Sfast) and the rate with which carbon is transferred from this stock to the atmosphere (tsfast_A).

            Returns:
                None.
        """

        self.Sfast["current"]=self.Sfast["previous"]+self.tsdead_Sfast*self.Sdead["previous"]-self.tsfast_A*self.Sfast["previous"]

    def calculate_NEE(self):
        """
        Return the Net Ecosystem Exchange (NEE) for the previous timestep.

        This is all the carbon that was absorbed by living trees (Cgpp) minus
        what was emitted to the atmosphere by the living trees (Cr) and
        the decomposition of the dead wood and soil stocks (tsdead_A+tsslow_A+tsfast_A).

        Returns:
            None.

        """

        #self.smort_log.append(self.Smort)

        self.calculate_Sdead()
        self.calculate_Sslow()
        self.calculate_Sfast()


        NEE= self.Cgpp-self.Cr-(self.tsdead_A*self.Sdead["previous"])-(self.tsslow_A*self.Sslow["previous"])-(self.tsfast_A*self.Sfast["previous"])

        self.Sdead["previous"]=self.Sdead["current"]
        self.Sslow["previous"]=self.Sslow["current"]
        self.Sfast["previous"]=self.Sfast["current"]



        return NEE
