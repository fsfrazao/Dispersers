from py_ibm import *
#from plot_trees import plot_tree_pos
import pdb
from tables import *
import json

class DispersersGrid(Rectangular_Grid):
    def __init__(self,x_max,y_max,patch_area,parent_tree_class,ndim=1,dim_names=["Resistance"], ):
        super().__init__(x_max,y_max,ndim=ndim,dim_names=dim_names)
        self.patch_area=patch_area
        self.total_area=self.x_max * self.y_max *self.patch_area/10000
        self.trees_per_patch={(x,y):[] for x in range(0,self.x_max) for y in range(0,self.y_max)}
        self.parent_tree_class=parent_tree_class
        self.FTs=[k for k in self.parent_tree_class.DERIVED_TREES.keys()]
        self.seedbank={k:[] for k in self.FTs}

    def update(self):
        """ Execute the methods that calculate values for each dimension of the grid.

        Returns:
            None.

        """
        empty_patches={p:0 if len(t)>0 else 1 for p,t in self.trees_per_patch.items()}
        for k,v in empty_patches.items():
            self.surface[self.dim_names["Resistance"]][k]=v

    def trees_in_patches(self,patches):
        """ Make a list of tree ids and their respective positions.

        Args:
            patch (list): a list of tuples representing the patch coordinates '(x,y)'.

        Returns:
            id_pos (list):a list tuples with tree ids and positions '(id,(x_pos,y_pos))'

        """


        #one list of ids for each patch
        ids=[self.trees_per_patch[p] for p in patches]
        #one single list of tuples with (id, pos)
        id_pos=[(tree,self.parent_tree_class.Instances.get(tree).position) for patch in ids for tree in patch]

        return id_pos
