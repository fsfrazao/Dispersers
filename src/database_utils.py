from tables import *

def open_database(file_name, simulation):
    h5_database=open_file(file_name,mode="a")

    sim_node="/sim_{0}".format(simulation)
    try:
        h5_database.get_node(sim_node)
    except NoSuchNodeError:
        sim=h5_database.create_group("/","sim_{0}".format(simulation),"Simulation {0}".format(simulation))

    return h5_database


def create_tables(h5_database,simulation=1):

    class Individuals(IsDescription):
        time_step=StringCol(4)
        day=StringCol(4)
        month=StringCol(4)
        year=StringCol(4)
        ind_id=Int32Col()
        initial_x=Float32Col()
        initial_y=Float32Col()
        final_x=Float32Col()
        final_y=Float32Col()
        energy=Float32Col()

    class Seeds(IsDescription):
        time_step=StringCol(4)
        day=StringCol(4)
        month=StringCol(4)
        year=StringCol(4)
        tree_id=Int32Col()
        initial_x=Float32Col()
        initial_y=Float32Col()
        final_x=Float32Col()
        final_y=Float32Col()
        dispersal_type=StringCol(24)


    class Populations(IsDescription):
        time_step=StringCol(4)
        day=StringCol(4)
        month=StringCol(4)
        year=StringCol(4)
        PrimaryDisperser=Int32Col()



    dispersers_lvl=h5_database.create_group("/sim_{0}".format(simulation),"dispersers","Dispersers Outputs for Simulation {0}".format(simulation))
    sys_lvl=h5_database.create_group("/sim_{0}/dispersers".format(simulation),"sys_lvl","System Level Observations for Simulation {0}".format(simulation))

    Pop_table=h5_database.create_table(sys_lvl,"Populations",Populations,"Disperser Population Sizes")
    
    Seeds_table=h5_database.create_table(sys_lvl,"Seeds",Seeds,"Seed Dispersal Data")

    ind_lvl=h5_database.create_group("/sim_{0}/dispersers".format(simulation),"ind_lvl","Individual Level Observations for Simulation {0}".format(simulation))

    Ind_table=h5_database.create_table(ind_lvl,"Ind",Individuals,"Individual Level Data")
