from scipy.spatial import ConvexHull
from pandas import read_hdf
import numpy as np

# points=np.random.rand(1000,4)*100
# df=pd.DataFrame(points, columns=["initial_x","initial_y", "final_x", "final_y"])
# df=df.apply(func=lambda x: round(x,2))
# ids=pd.DataFrame(np.random.choice(range(10), size=1000, replace=True),columns=["id"])
# dates=pd.DataFrame(np.random.choice(["1","2","3"], size=1000, replace=True), columns=["date"])
# df=df.join([ids,dates])
#
# df["path_length"]=df.apply(calc_dist,axis=1)
# path_length_by_ind=df.groupby("id").aggregate(np.mean)
# avg_path_length=np.mean(path_length_by_ind.path_length)

def calc_dist(row):
     return  euclid_dist(row.initial_x,row.initial_y,row.final_x,row.final_y)

def area_final_positions(data):
    try:
        hull=ConvexHull(data[["final_x","final_y"]])
        return hull.area
    except:
        return float('nan')

def euclid_dist(x0,y0,x1,y1):
    return np.sqrt((x1-x0)**2+(y1-y0)**2)


def hdf_table_to_df(database,table):
    read_hdf(database,)



def avg_distance(data):
    data["distance"]=data.apply(calc_dist,axis=1)
    #path_length_by_ind=df.groupby("id").aggregate(np.mean)
    avg_distance=np.mean(data.distance)
    avg_distance=round(float(avg_distance),2)

    return avg_distance



def avg_area(data, ind_key):
    ind_group=data.groupby(ind_key)
    areas=ind_group.apply(area_final_positions)
    avg_area=np.mean(areas)
    return avg_area



def home_range(database, sim_number):
    #TODO: add START and END parameters
    data=read_hdf(database, "/sim_{0}/dispersers/ind_lvl/Ind".format(sim_number))
    home_range=avg_area(data,"ind_id")

    return home_range


def avg_path_length(database,sim_number):
    #TODO: add START and END parameters
    data=read_hdf(database, "/sim_{0}/dispersers/ind_lvl/Ind".format(sim_number))
    avg_path_length=avg_distance(data)

    return avg_path_length


def seed_shadow(database,sim_number):
    #TODO: add START and END parameters
    data=read_hdf(database, "/sim_{0}/dispersers/sys_lvl/Seeds".format(sim_number))
    seed_shadow=avg_area(data,"tree_id")

    return seed_shadow


def avg_dispersal_distance(database,sim_number):
    #TODO: add START and END parameters
    data=read_hdf(database, "/sim_{0}/dispersers/sys_lvl/Seeds".format(sim_number))
    avg_dispersal=avg_distance(data)

    return avg_dispersal
