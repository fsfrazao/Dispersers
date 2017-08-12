from py_ibm import *
from plot_trees import plot_tree_pos
import pdb
from tables import *
import json


class Dispersers_World(Tree_World):

    def __init__(self,topology, parent_tree_class):
        super().__init__(topology)
        self.seeds_db=None
        self.seeds_cursor=None
        self.parent_tree_class=parent_tree_class

    def open_seed_database(self, path="./"):
        self.seeds_db=connect(path+"seed_db.sql")
        self.seeds_cursor=self.seeds_db.cursor()
        self.seeds_cursor.execute(''' CREATE TABLE IF NOT EXISTS seeds(tree_id INTEGER, initial_x FLOAT, initial_y FLOAT, final_x FLOAT, final_y FLOAT, dispersal_type TEXT)''')
        self.seeds_db.commit()

    def add_seed_entry(self,tree_id,initial_x,
        initial_y,final_x,final_y,dispersal_type):

        self.seeds_cursor.execute(''' INSERT INTO seeds(tree_id,initial_x,initial_y,final_x,final_y,dispersal_type) VALUES(?,?,?,?,?,?)''',(tree_id,initial_x,
            initial_y,final_x,final_y,dispersal_type))
        #self.seeds_db.commit()

    def close_seeds_db(self):
        self.seeds_db.commit()
        self.seeds_db.close()

    def increase_dispersers_population(self):
        n=PrimaryDisperser.CalculatePopGrowth()
        self.create_dispersers(n)

    def create_dispersers(self,n):
        trees_for_dispersers=np.random.choice(list(self.parent_tree_class.Instances.keys()),size=n,replace=True)

        for i in range(n):
            tree_id=trees_for_dispersers[i]
            pos=Tree.Instances.get(tree_id).position
            PrimaryDisperser(position=pos,world=self,gut_passage_time=2, max_fruits=5, current_tree_id=tree_id)
