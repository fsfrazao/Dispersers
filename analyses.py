from scipy.spatial import ConvexHull
import pandas as pd
import numpy as np

points=np.random.rand(1000,4)*100
df=pd.DataFrame(points, columns=["initial_x","initial_y", "final_x", "final_y"])
df=df.apply(func=lambda x: round(x,2))
ids=pd.DataFrame(np.random.choice(range(10), size=1000, replace=True),columns=["id"])
dates=pd.DataFrame(np.random.choice(["1","2","3"], size=1000, replace=True), columns=["date"])
df=df.join([ids,dates])

df["path_length"]=df.apply(calc_dist,axis=1)
path_length_by_ind=df.groupby("id").aggregate(np.mean)
avg_path_length=np.mean(path_length_by_ind.path_length)

def calc_dist(row):
     return  euclid_dist(row.initial_x,row.initial_y,row.final_x,row.final_y)

def area_final_positions(data):
    return ConvexHull(data[["final_x","final_y"]]).area




def euclid_dist(x0,y0,x1,y1):
    return np.sqrt((x1-x0)**2+(y1-y0)**2)


def avg_distance(data):
    data["path_length"]=data.apply(calc_dist,axis=1)
    #path_length_by_ind=df.groupby("id").aggregate(np.mean)
    avg_path_length=np.mean(data.path_length)
    avg_path_length=round(float(avg_path_length))

    return avg_path_length



def avg_area(data, ind_key)):
    ind_group=data.groupby(ind_key)
    areas=ind_group.apply(area_final_positions)
    avg_area=np.mean(areas)
    return avg_area



def home_range(database, sim, from=None, to=None):
    pass

def avg_path_length():
    pass


def dispersal_distance():
    pass

def seed_shadow():
    pass
