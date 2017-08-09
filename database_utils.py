from tables import *

def create_HDF_database(self,database_name):

    class Individuals(IsDescription):
        time_stamp=StringCol()
        ind_id=Int32Col()
        pos_x=Float32Col()
        pos_y=Float32Col()
        energy=Float32Col()

    class Seeds(IsDescription):
        deposition_time=StringCol()
        tree_id=Int32Col()
        initial_x=Float32Col()
        initial_y=Float32Col()
        final_x=Float32Col()
        final_y=Float32Col()
        dispersal_type=StringCol()


    class Population(IsDescription):
        time_stamp=StringCol()
        PrimaryDisperser=Int32Col()



    h5file=open_file(database_name,mode="w",title="Tree Model outputs")

    sim=h5file.create_group("/","sim_1","Simulation 1")
    sys_lvl=h5file.create_group("/sim_1","sys_lvl","System Level Observations for simulation 1")

    NEE_table=h5file.create_table(sys_lvl,"NEE",NEE,"Net Ecosystem Exchange")
    Pop_table=h5file.create_table(sys_lvl,"Pop",Populations,"Population Sizes")
    Stocks_table=h5file.create_table(sys_lvl,"Stocks",CarbonStocks,"Carbon Stocks")
    #Emissions_table=h5file.create_table(sys_lvl,"Emissions",CarbonEmissions,"Emissions from soil stocks, living trees and dead wood")

    ind_lvl=h5file.create_group("/sim_1","ind_lvl","Individual Level Observations for Simulation 1")

    Ind_table=h5file.create_table(ind_lvl,"Ind",Individuals,"Individual Level Data")

    h5file.close()
