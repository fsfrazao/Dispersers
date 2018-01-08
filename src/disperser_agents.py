from py_ibm import *
from trees_ibm import *
from math import floor
import numpy as np



class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class TreeError(Error):
    """Exception raised for errors regarding trees.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message



class PrimaryDisperser(Agent):

    GrowthCoefficient=0.1
    EnergyPerFruit=8
    CostPerDisperser=60
    AvgEnergy=880

    DERIVED_DISPERSERS={}

    @classmethod
    def DisperserFactory(cls,new_cls_name,new_parameters):

        def new_init(self,position,world,id=None, **kwargs):



            super().__init__(position,world,id=id)


            self.patch=(floor(self.position[0]), floor(self.position[1]))
            self.seeds=[]
            self.track=[]
            self.energy_track=[]
            self.activity_track=[]
            self.time_on_this_tree=0
            self.previous_position=self.position

            self.gut_passage_time=new_parameters['gut_passage_time']
            self.max_fruits=new_parameters['max_fruits']
            self.current_tree_id=new_parameters['current_tree_id']
            self.energy=new_parameters['energy']
            self.min_rep_energy=new_parameters['min_rep_energy']
            self.energy_level_1=new_parameters['energy_level_1']
            self.energy_level_2=new_parameters['energy_level_2']
            self.tree_time=new_parameters['tree_time']
            self.feeding=new_parameters['feeding']
            self.traveling=new_parameters['traveling']
            self.resting=new_parameters['resting']
            self.radius=new_parameters['radius']

            PrimaryDisperser.IncrementID()


            if id == None:
                #self.id=Tree.ID
                self.Indices.append(self.id)
                cls.Instances[self.id]=self
                cls.IncrementID()
            else:
                self.id=id
                self.Indices.append(self.id)
                cls.Instances[self.id]=self

        new_attr_met={"Type":new_cls_name,
                    'EnergyPerFruit':new_parameters['EnergyPerFruit'],
                    'CostPerDisperser':new_parameters['CostPerDisperser'],
                    'AvgEnergy':new_parameters['AvgEnergy'],
                    'Indices':[],
                    '__init__':new_init}


        new_cls=type(new_cls_name, (cls,), new_attr_met)

        cls.DERIVED_DISPERSERS[new_cls_name]=new_cls

        return new_cls

    @classmethod
    def TotalEnergy(cls):
        tEnergy=0
        count=0
        for i in cls.Instances.values():
            if i.energy >= i.min_rep_energy:
                tEnergy+=i.energy-i.min_rep_energy
                count+=1

        return (count,tEnergy)

    @classmethod
    def PopulationSize(cls):
        return len(cls.Instances)

    @classmethod
    def CalculatePopGrowth(cls,world):
        N_r, Et= cls.TotalEnergy()
        N_0=cls.PopulationSize()
        El=world.parent_tree_class.CountFruits()*cls.EnergyPerFruit
        K=El/cls.AvgEnergy
        alpha=1/cls.CostPerDisperser
        growth=alpha*Et*((K-N_0)/K) if K>0 else 0
        #growth_rate=cls.TotalEnergy()*cls.GrowthCoefficient
        #growth=N_r*growth_rate
        return int(round(growth))

    def __init__(self,position,world, gut_passage_time,max_fruits, current_tree_id=None):
        super().__init__(position,world)
        self.gut_passage_time=gut_passage_time
        self.max_fruits=max_fruits
        self.seeds=[]
        self.patch=(floor(self.position[0]), floor(self.position[1]))
        self.energy=60
        self.min_rep_energy=70
        self.energy_level_1=80
        self.energy_level_2=150
        self.tree_time=5
        self.time_on_this_tree=0
        self.feeding=8
        self.traveling=-1.6
        self.resting=-0.7
        self.previous_position=self.position
        self.current_tree_id=current_tree_id
        self.rest_prob=0.5
        self.radius=50
        self.track=[]
        self.energy_track=[]
        self.activity_track=[]
        PrimaryDisperser.IncrementID()


    def increase_tree_timer(self):
        self.time_on_this_tree+=1

    def reset_tree_timer(self):
        self.time_on_this_tree=0

    def update_position(self,pos):
        new_position=self.world.topology.in_bounds(pos)
        self.previous_position=self.position
        self.position=new_position
        self.patch=(floor(self.position[0]), floor(self.position[1]))

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
        #verified_position=self.world.topology.in_bounds(new_position)
        #self.update_position(verified_position)
        self.update_position(new_position)
        #self.track.append(self.position)

    def move(self,pos):

        self.update_position(pos)
        self.energy+=self.traveling
        #self.track.append(self.position)

    def eat(self,tree_id):
        tree=self.world.parent_tree_class.Instances.get(tree_id)
        tree_pos=tree.position
        tree_type=tree.Ftype
        self.add_seed(tree_pos=tree_pos,tree_type=tree_type,tree_id=tree_id)
        self.energy+=self.feeding

        tree.decrement_fruit_stock()

    def digest_seeds(self):
        self.seeds=sorted(self.seeds, key=lambda k: k['time'])
        for i,s in enumerate(self.seeds):
            if s["time"]<(self.world.step-self.gut_passage_time):
                self.defecate(self.seeds.pop(i))

    def defecate(self,seed):
        #seed_pos=(self.position[0]+random()*0.5, self.position[1]+random()*0.5 )
        x0,y0=self.previous_position
        x1,y1=self.position
        seed_pos=self.world.topology.rand_point_along_line(x0=round(x0,2),y0=round(y0,2),x1=round(x1,2),y1=round(y1,2),offset=2)
        seed_pos=self.world.topology.in_bounds(seed_pos)
        self.world.topology.seedbank[seed['tree_type']].append(seed_pos)
        self.world.add_seed_entry(tree_id=seed["tree_id"],initial_x=seed["tree_pos"][0],
            initial_y=seed["tree_pos"][1],final_x=seed_pos[0],final_y=seed_pos[1],dispersal_type="primary")


    def add_seed(self, tree_pos,tree_type,tree_id):
        """ Adds a seed to the disperser's gut"""
        self.seeds.append({"time":self.world.step,"tree_type":tree_type,"tree_id":tree_id,"tree_pos":tree_pos})

    def trees_within_radius(self,radius):
        cell_length=self.world.topology.patch_area**0.5

        neighborhood=self.world.topology.hood(cell=self.patch,remove_center=False,radius=round(radius/cell_length))
        trees_in_neighborhood=self.world.topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if self.world.topology.point_within_circle(self.position[0],self.position[1],radius,p[0],p[1],self.world.topology.scaled_distance,scale=cell_length)]

        if len(trees_in_circle)>0:
            return trees_in_circle
        else:
            raise TreeError(None,"no tree within radius")

    def feeding_tree_value(self,tree):
        resistance_penalty=1.2
        distance=self.world.topology.euclidean_distance(self.position[0],self.position[1],tree.position[0],tree.position[1])

        gap=self.world.topology.line_between_points(self.patch[0],self.patch[1],tree.patch[0],tree.patch[1])
        x=[g[0] for g in gap]
        y=[g[1] for g in gap]
        resistance=sum(self.world.topology.surface[0][x,y])
        resistance=resistance*resistance_penalty
        if resistance>0:
            distance+=resistance
        value=tree.available_fruits/distance if distance !=0 else 0


        return value

    def non_feeding_tree_value(self,tree):
        resistance_penalty=0.3
        distance=self.world.topology.euclidean_distance(self.position[0],self.position[1],tree.position[0],tree.position[1])

        gap=self.world.topology.line_between_points(self.patch[0],self.patch[1],tree.patch[0],tree.patch[1])
        x=[g[0] for g in gap]
        y=[g[1] for g in gap]
        resistance=sum(self.world.topology.surface[0][x,y])
        resistance=resistance*resistance_penalty
        if resistance>0:
            distance+=resistance

        value=1/distance if distance !=0 else 0


        return value



    #TODO:
    #What if there are no trees around?
    def choose_feeding_tree(self,trees):
        if len(trees)>0:
            tree_values=[]
            for tree in trees:
                tree_value=self.feeding_tree_value(self.world.parent_tree_class.Instances.get(tree))
                tree_values.append({'tree_id':tree,'tree_value':tree_value})

            sorted_by_value=sorted(tree_values, key=lambda k: k['tree_value'],reverse=True)

            return(sorted_by_value[0]['tree_id'])
        else:
            raise TreeError(None,"no feeding tree")

    def choose_non_feeding_tree(self,trees):
        if len(trees)>0:
            tree_values=[]
            for tree in trees:
                tree_value=self.non_feeding_tree_value(self.world.parent_tree_class.Instances.get(tree))
                tree_values.append({'tree_id':tree,'tree_value':tree_value})

            sorted_by_value=sorted(tree_values, key=lambda k: k['tree_value'],reverse=True)

            return(sorted_by_value[0]['tree_id'])
        else:
            raise TreeError(None,"no feeding tree")


    def active_or_not(self):
        if np.random.rand() < self.rest_prob:
            self.active=False
        else:
            self.active=True

    def on_feeding_tree(self):

        if self.world.parent_tree_class.Instances.get(self.current_tree_id).available_fruits>0:
            return True
        else:
            return False


    def random_activity(self):
        p=(0.5,0.5)
        activity=np.random.choice(('resting','moving'),p=p)
        if activity=='resting':
            self.activity_resting()
        else:
            self.activity_roam()


    def activity_feeding(self):
        self.activity_track.append("feed")
        if self.on_feeding_tree():
            self.eat(self.current_tree_id)

    def activity_go_to_feeding_tree(self):
        try:
            trees=self.trees_within_radius(self.radius)
        except TreeError:
            print("There are no other trees around")
        trees.remove(self.current_tree_id)

        try:
            chosen_tree=self.world.parent_tree_class.Instances.get(self.choose_feeding_tree(trees))
            self.move(pos=chosen_tree.position)
            self.current_tree_id=chosen_tree.id
        except:
            print("no other trees")


        #What if the only tree around is the current tree?

    def activity_resting(self):
        self.activity_track.append("rest")
        self.energy+=self.resting
        self.increase_tree_timer()


    def activity_roam(self):
        self.activity_track.append("roam")
        trees=self.trees_within_radius(self.radius)

        if len(trees) >=1:
            chosen_tree=self.world.parent_tree_class.Instances.get(self.choose_non_feeding_tree(trees))#self.world.parent_tree_class.Instances.get(np.random.choice(trees))
            self.move(pos=chosen_tree.position)
            self.current_tree_id=chosen_tree.id
            self.reset_tree_timer()
            self.increase_tree_timer()
        #What if the only tree around is the current tree?

    def low_energy_action(self):
        if self.on_feeding_tree():
            if self.time_on_this_tree< self.tree_time:
                self.activity_feeding()
                self.increase_tree_timer()
            else:
                self.activity_go_to_feeding_tree()
                self.reset_tree_timer()
                self.activity_feeding()
                self.increase_tree_timer()
        else:
            self.activity_go_to_feeding_tree()
            self.reset_tree_timer()
            self.activity_feeding()
            self.increase_tree_timer()

    def moderate_energy_action(self):
        self.random_activity()

    def high_energy_action(self):
        self.activity_resting()

    def schedule(self):
        # if self.energy < self.energy_level_1:
        #     self.low_energy_action()
        # else:
        #     if self.from_feeding():
        #         if self.energy>self.energy_level_2:
        #             self.high_energy_level()
        #         else:
        #             self.low_energy_action()
        #     else:
        #         self.energy+=self.resting
        self.track.append(self.position)
        if self.energy<self.energy_level_2:
            trees=self.trees_within_radius(50)
            chosen_tree=self.world.parent_tree_class.Instances.get(self.choose_feeding_tree(trees))
            if chosen_tree != None:
                self.move(pos=chosen_tree.position)
                self.current_tree_id=chosen_tree.id
                self.eat(tree=chosen_tree)
                #print(d.energy, ">>>>>>")
            else:
                self.activity_resting()
                #print(d.energy, "------|")
        else:
            self.activity_resting()
            #print(d.energy, "------")
        self.digest_seeds()

    def remove_me(self):
        del PrimaryDisperser.Instances[self.id]
    
    def current_tree_exists(self):
        '''Checks if the current tree still exists'''
        return self.world.parent_tree_class.Instances.get(self.current_tree_id) is not None

    def schedule2(self):
        self.energy_track.append(self.energy)
        self.track.append(self.position)
        
        #if the current tree does not exist anymore,
        # go to a new tree
        if not self.current_tree_exists():
            self.activity_roam()
            self.energy-=self.traveling # add the spent energy again
                    

        if self.energy<=0:
            self.remove_me()
        elif self.energy<=self.energy_level_1:
            self.low_energy_action()
        elif self.energy<=self.energy_level_2:
            self.moderate_energy_action()
        else:
            self.high_energy_action()


        self.digest_seeds()
