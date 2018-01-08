from py_ibm import *
#from plot_trees import plot_tree_pos
#import pdb
from tables import *
import json
from numpy import floor


class SimpleTree(Agent):
    ID = 0
    Instances = {}
    DERIVED_TREES = {}

    @classmethod
    def TreeFactory(cls, new_cls_name, new_parameters):

        def new_init(self, position, world, id=None, **kwargs):

            super().__init__(position, world, id=id)

            self.patch = (int(floor(self.position[0])), int(
                floor(self.position[1])))
            self.available_fruits = 0
            self.world.topology.trees_per_patch[self.patch].append(self.id)

            if id == None:
                # self.id=Tree.ID
                self.Indices.append(self.id)
                cls.Instances[self.id] = self
                cls.IncrementID()
            else:
                self.id = id
                self.Indices.append(self.id)
                cls.Instances[self.id] = self

        new_attr_met = {  # 'Nseed':new_parameters['Nseed'],
            'Indices': [],
            'Ftype': new_cls_name,
            '__init__': new_init}

        new_cls = type(new_cls_name, (cls,), new_attr_met)

        cls.DERIVED_TREES[new_cls_name] = new_cls

        return new_cls

    @classmethod
    def CountFruits(cls):
        total_fruits = 0
        for t in cls.Instances.values():
            total_fruits += t.available_fruits
        return total_fruits

    def __repr__(self):
        return "Tree id: {0}".format(self.id)

    def __init__(self, position, world, id=None):

        super().__init__(position, world, id=id)
        # self.AddInstance(self.id)
        #Tree.Instances[Ftype+' '+str(self.id)]=self

        self.patch = (int(floor(self.position[0])), int(
            floor(self.position[1])))
        self.availbale_fruits = 0
        self.world.topology.trees_per_patch[self.patch].append(self.id)

    def decrement_fruit_stock(self, n=1):
        self.available_fruits -= n

    def disperse_seeds(self, fraction=0.5, a=0.3):
        n_seeds = int(round(self.available_fruits * fraction))
        self.decrement_fruit_stock(n=n_seeds)
        seed_distances = np.random.power(a=a, size=n_seeds)
        for d in seed_distances:
            seed_position = self.random_seed_position(distance=d)
            self.deposit_seed(seed_position)

    def random_seed_position(self, distance):
        x0, y0 = self.position
        angle = np.random.randint(360)
        new_position = ((x0 + cos(radians(angle)) * distance,
                         y0 + sin(radians(angle)) * distance))
        verified_position = self.world.topology.in_bounds(new_position)
        return verified_position

    def deposit_seed(self, seed_position):
        self.world.topology.seedbank[self.Ftype].append(seed_position)
        self.world.add_seed_entry(tree_id=self.id, initial_x=self.position[0],
                                  initial_y=self.position[1], final_x=seed_position[0], final_y=seed_position[1], dispersal_type="non-zoochoric")
