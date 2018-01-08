from py_ibm import *
from trees_ibm.world import Tree_World
from dispersers import database_utils #import dispersers.database_utils
from dispersers.disperser_agents import PrimaryDisperser
#from plot_trees import plot_tree_pos
from tables import *
import json
import numpy as np



class Dispersers_World(World):

    def __init__(self,topology, parent_tree_class):
        super().__init__(topology)
        self.db=None
        self.seed_table=None
        self.ind_table=None
        self.sim_number=1
        self.parent_tree_class=parent_tree_class

    def setup_database(self, filename, sim_number):
        self.db=database_utils.open_database(filename,sim_number)
        try:
            self.db.get_node("/sim_{0}/dispersers/sys_lvl/Seeds".format(sim_number))
        except NoSuchNodeError:
            database_utils.create_tables(self.db,sim_number)

        self.seed_table=self.db.get_node("/sim_{0}/dispersers/sys_lvl/Seeds".format(sim_number))
        self.ind_table

    def add_seed_entry(self,tree_id,initial_x,
        initial_y,final_x,final_y,dispersal_type):
        seed_table=self.db.get_node("/sim_{0}/dispersers/sys_lvl/Seeds".format(self.sim_number))
        seed_r=seed_table.row

        seed_r["time_step"]=self.step
        seed_r["day"]=self.day
        seed_r["month"]=self.month
        seed_r["year"]=self.year
        seed_r["tree_id"]=tree_id
        seed_r["initial_x"]=initial_x
        seed_r["initial_y"]=initial_y
        seed_r["final_x"]=final_x
        seed_r["final_y"]=final_y
        seed_r["dispersal_type"]=dispersal_type

        seed_r.append()
        seed_table.flush()

    def add_pop_entry(self):
        pop_table=self.db.get_node("/sim_{0}/dispersers/sys_lvl/Populations".format(self.sim_number))
        pop_r=pop_table.row
        pop_r["time_step"]=self.step
        pop_r["day"]=self.day
        pop_r["month"]=self.month
        pop_r["year"]=self.year
        pop_r['PrimaryDisperser']=PrimaryDisperser.PopulationSize()
        
        pop_r.append()
        pop_table.flush()
        
    def add_ind_entry(self,ind_id,initial_x,
        initial_y,final_x,final_y,energy):
        ind_table=self.db.get_node("/sim_{0}/dispersers/ind_lvl/Ind".format(self.sim_number))
        ind_r=ind_table.row

        ind_r["time_step"]=self.step
        ind_r["day"]=self.day
        ind_r["month"]=self.month
        ind_r["year"]=self.year
        ind_r["ind_id"]=ind_id
        ind_r["initial_x"]=initial_x
        ind_r["initial_y"]=initial_y
        ind_r["final_x"]=final_x
        ind_r["final_y"]=final_y
        ind_r["energy"]=energy

        ind_r.append()
        ind_table.flush()


    def run_simulation(self,n):
        #self.db=database_utils.create_database(database_name,sim_number)

        for i in range(n):
            dispersers_ids=list(PrimaryDisperser.Instances.keys())
            for i in dispersers_ids:
                d=PrimaryDisperser.Instances.get(i)
                d.schedule2()
                self.add_ind_entry(ind_id=d.id,
                initial_x=d.previous_position[0],
                initial_y=d.previous_position[1],
                final_x=d.position[0],
                final_y=d.position[1],
                energy=d.energy)
            self.increment_time()
        self.db.close()




    def close_seeds_db(self):
        self.seeds_db.commit()
        self.seeds_db.close()

    def increase_dispersers_population(self):
        n=PrimaryDisperser.CalculatePopGrowth(world=self)
        self.create_dispersers(n)

    def create_dispersers(self,n, disperser_class=PrimaryDisperser):
        trees_for_dispersers=np.random.choice(list(self.parent_tree_class.Instances.keys()),size=n,replace=True)

        for i in range(n):
            tree_id=trees_for_dispersers[i]
            pos=self.parent_tree_class.Instances.get(tree_id).position
            disperser_class(position=pos,world=self,gut_passage_time=2, max_fruits=5, current_tree_id=tree_id)


    def create_trees(self,tree_class,n,pos=None, ids=None, **kwargs):
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
                tree_class(position=pos[i],**kwargs)

        else:
            for i in range(n):
                tree_class(position=pos[i], id=ids[i],**kwargs)




    def seedbank_from_file(self, input_file):
        with open(input_file,'r') as json_file:
            self.topology.seedbank=json.load(json_file)



    def seedbank_to_file(self,output_file):
        with open(output_file,'w') as json_file:
            json.dump(self.topology.seedbank,json_file)


    def model_status_from_file(self, input_file):
        with open(input_file,'r') as json_file:
            trees=json.load(json_file)

        for ft in self.parent_tree_class.DERIVED_TREES.keys():
            positions=[]
            ids=[]
            for tree in trees[ft]:
                positions.append(tuple(tree['position']))
                ids.append(tree['id'])

            self.create_trees(self.parent_tree_class.DERIVED_TREES[ft],len(trees[ft]),pos=positions,ids=ids, world=self)

        self.topology.update()
        # for t in Tree.Instances.values():
        #     t.update()
        #     t.AGB=t.AGB_from_DBH(t.DBH/100)
        self.parent_tree_class.ID=max(self.parent_tree_class.Instances.keys())+1


    def model_status_to_file(self,output_file):
        trees_per_type={k:[] for k in self.parent_tree_class.DERIVED_TREES.keys()}
        for t_id,t in self.parent_tree_class.Instances.items():
            trees_per_type[t.Ftype].append({"id":t_id, "position":t.position})

        with open(output_file,'w') as json_file:
            json.dump(trees_per_type,json_file)
