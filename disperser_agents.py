from py_ibm import *
from trees_ibm import *


class PrimaryDisperser(Agent):

    def __init__(self,position,world, gut_passage_time,max_fruits):
        super().__init__(position,world)
        self.gut_passage_time=gut_passage_time
        self.max_fruits=max_fruits
        self.seeds=[]



    def move_forward(self,distance):
        """ Move beetle forward in a straight line.

        Args:
            distance (float): distance of movement in grid cells.

        Returns:
                None.
        """

        x0,y0=self.position
        angle=self.orientation
        new_position=((x0+cos(radians(angle))*distance,y0+sin(radians(angle))*distance))
        verified_position=self.world.topology.in_bounds(new_position)
        self.update_position(verified_position)
        self.track.append(self.position)

    def move(self,pos):
        self.update_position(pos)

    def eat(self,tree):
        tree_pos=tree.position
        tree_type=tree.Ftype
        self.add_seed(tree_pos=tree_pos,tree_type=tree_type)
        tree.decrement_fruit_stock()

    def digest_seed(self):
        sorted_by_time=sorted(self.seeds, key=lambda k: k['time'])
        for i in range(len(sorted_by_time)):
            if sorted_by_time[i]["time"]>(self.world.step-self.gut_passage_time):
                self.defecate(sorted_by_time.pop(i))

    def defecate(self,seed):
        pass


    def add_seed(self, tree_pos,tree_type):
        self.seeds.append({"time":self.world.step,"tree_type":tree_type,"tree_pos":tree_pos})

    def trees_within_radius(self,radius):

        neighborhood=self.world.topology.hood(cell=self.patch,remove_center=False)
        trees_in_neighborhood=self.world.topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if self.world.topology.point_within_circle(self.position[0],self.position[1],radius,p[0],p[1],self.world.topology.scaled_distance,scale=self.world.topology.patch_area**0.5)]

        return trees_in_circle

    def feeding_tree_value(self,tree):
        distance=self.world.topology.euclidean_distance(self.position[0],self.position[1],tree.position[0],tree.position[1])

        value=distance/tree.available_fruits
        return value

    def choose_feeding_tree(self,trees):
        tree_values=[]
        for tree in trees:
            tree_value=self.feeding_tree_value(Tree.Instances.get(tree))
            tree_values.append({'tree_id':tree,'tree_value':tree_value})

        sorted_by_value=sorted(tree_values, key=lambda k: k['tree_value'],reverse=True)

        return(sorted_by_value[0]['tree_id'])
