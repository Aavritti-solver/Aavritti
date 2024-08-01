'''
settings.py
===============
created : Vighnesh JR (27/7/2024)
contains class 
- Settings
'''
class Settings:
    '''
    Settings
    ========
    created: Vighnesh JR on 27/7/2024
    ----------------------------------
    student at IIT Bombay
    ---------------------
    Class that stores the following 
    1. directory information
    2. Constants (dictionary format) 
    3. Computational flag (boolean)
    4. Boundary condition details
    5. Element specifications
    '''
    def __init__(self,
                 mesh_folder      = None,
                 mesh_name        = None,
                 baseflow_folder  = None,
                 baseflow_name    = None,
                 result_folder    = None,
                 result_name      = None,
                 constants        = None,
                 values           = None,
                 boundary         = None,

                 # Default case for optimal computation
                 computation_flag = False,
                 is_high_order    = False 
                 ):

        # Directories
        self.mesh_folder        = mesh_folder
        self.mesh_name          = mesh_name
        self.baseflow_folder    = baseflow_folder
        self.baseflow_name      = baseflow_name
        self.result_folder      = result_folder
        self.result_name        = result_name

        # Replace with dictonary
        self.constants          = constants
        self.values             = values
        self.boundary           = boundary

        # Boolean flags
        self.computaion_flag    = computation_flag
        self.is_high_order      = is_high_order        
    