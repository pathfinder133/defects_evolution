# #Superkomórki - komórki które leżą poza siatką i są używane do wyznaczania scieżeki defektu (szczeliny)
# #jako droga początkowa i końcowa sćieżki
# #
# #                                        SUPERCELL_F
# #                         #################################
# #                       #                               # #
# #                     #          SUPERCELL_A          #   #
# #                  ##################################     #
# #                  #                                #     #
# #                  #                                #     #
# #  SUPERCELL_D     #                                #  SUPERCELL_B
# #                  #  SUPERCELL_E                   #     #
# #                  #                                #   #
# #                  #                                # #
# #                  ##################################
# #                               SUPERCELL_C
# #
# #
# # W przypadku  siatki w przestrzeni  SUPERCELL_E to kormórka przed ssześcianem, a SUPERCELL_F to komórka za sześcianem
# #


import json
from   dataclasses import dataclass, asdict
from   typing import Final, Self

ITERATIONS            = "ITERATIONS"

ENABLED               = "enabled"
DISABLED              = "disabled"

P1_PROBABILITY        = "P1_PROBABILITY"
P2_PROBABILITY        = "P2_PROBABILITY"

GRID_SIZE             = "GRID_SIZE"
GRID_SIZE_PIXELS      = "GRID_SIZE_PIXELS"
Z_AXIS                = "Z_AXIS"

HEAL_CYKLES           = "HEAL_CYKLES"
FPS                   = "FPS"

P1_COLOR              = "P1_COLOR"
P2_COLOR              = "P2_COLOR"
D_COLOR               = "D_COLOR"
DP_COLOR              = "DP_COLOR"
FP_COLOR              = "FP_COLOR"
PD_COLOR              = "PD_COLOR"

DEFECT                = "DEFECT"

SIM2D                 = "SIM2D"
SIM3D                 = "SIM3D"

SUPERCELL_A           = "SUPERCELL_A"
SUPERCELL_B           = "SUPERCELL_B"
SUPERCELL_C           = "SUPERCELL_C"
SUPERCELL_D           = "SUPERCELL_D"
SUPERCELL_E           = "SUPERCELL_E"
SUPERCELL_F           = "SUPERCELL_F"

DESTRUCTION_PATH      = "DESTRUCTION_PATH"
P1                    = "P1"
P2                    = "P2"
PERMANENT             = "PERMANENT"
FRACTURE_PATH         = "FRACTURE_PATH"
EDGE                  = "EDGE"

GRID_MARGINS          = "GRID_MARGINS"

DRYRUN                = "DRYRUN"
INCLUDE_CHANGED_CELLS = "INCLUDE_CHANGED_CELLS"
SIMULATION_CYKLES     = "SIMULATION_CYKLES"
#TYPE                  = "TYPE"
IMAGE_PATH            = "IMAGE_PATH"
BACKGROUND_COLOR      = "BACKGROUND_COLOR"



@dataclass
class Settings:
    """
    Holder for a various settings shared across all 
    modules. 

    It maintains settings, writing and recreating
    can be used as an object inection 
    """

    # How many simulation run one-by-one
    ITERATIONS: int

    # Draw probability for P1 and P2 cells
    P1_PROBABILITY :float
    P2_PROBABILITY :float

    # In-memory grid size in tiles
    GRID_SIZE :tuple[int,int]
    # Simulation_window size in pixels
    GRID_SIZE_PIXELS :tuple[int,int]
    # 3D in-memory cube representation, edge length
    # TODO: Rename it to something more corresponding like: EDGE3D_LENGTH
    #Z_AXIS : int

    # How many cykles will D cell have to heal
    HEAL_CYKLES :int
    # Simulation target FPS
    FPS :int
    #Color codes for coresponding cell types in visualisations
    P1_COLOR : tuple[int,int,int]
    P2_COLOR : tuple[int,int,int]
    D_COLOR  : tuple[int,int,int]
    DP_COLOR : tuple[int,int,int]
    FP_COLOR : tuple[int,int,int]
    PD_COLOR : tuple[int,int,int]
    BACKGROUND_COLOR: tuple[int, int, int]

    # Codes for supercells 
    SUPERCELL_A : int
    SUPERCELL_B : int
    SUPERCELL_C : int
    SUPERCELL_D : int
    SUPERCELL_E : int
    SUPERCELL_F : int

    # Number codes of partucular cell type
    DESTRUCTION_PATH : int
    P1 : int
    P2 : int
    DEFECT : int
    PERMANENT : int
    FRACTURE_PATH : int
    EDGE : int

    # Margins in px, used by  pygame.Surface.blit() method
    GRID_MARGINS : tuple[int,int]

    # Dryrun is without visualisation mode
    DRYRUN: int

    # Information for Data_collector() that it should include also
    # information from self.__changed if the simulation module
    #INCLUDE_CHANGED_CELLS: int
    # How many cykles will single simulation run (0 = inf)
    SIMULATION_CYKLES: int

    # Path to img file
    IMAGE_PATH: str


    @classmethod
    def from_defaults(cls, make_settings_file: bool = True) -> Self:
        """
        Default settings in the case, settings file is missing or
        something went wrong with the current settings file.
        """

        settings: Final[dict] = dict(

            P1_PROBABILITY        = 0.001,
            P2_PROBABILITY        = 0.002,
            P1                    = 1,
            P2                    = 10,
            DEFECT                = 1000,
            PERMANENT             = 10000,
            FRACTURE_PATH         = -200,
            DESTRUCTION_PATH      = -100,
            HEAL_CYKLES           = 5,
            EDGE                  = 100,
            P1_COLOR              = (21, 21, 21),
            P2_COLOR              = (50, 50, 50),
            D_COLOR               = (183, 28, 28),
            PD_COLOR              = (28, 255, 28),
            FP_COLOR              = (200, 200, 255),
            DP_COLOR              = (28, 28, 183),
            BACKGROUND_COLOR      = (0, 0, 0),
            SUPERCELL_A           = 1000010,
            SUPERCELL_B           = 1000020,
            SUPERCELL_C           = 1000030,
            SUPERCELL_D           = 1000040,
            SUPERCELL_E           = 1000050,
            SUPERCELL_F           = 1000060,
            GRID_SIZE             = (50,50),
            GRID_SIZE_PIXELS      = (1000, 800),
            #Z_AXIS                = 10,
            FPS                   = 10,
            GRID_MARGINS          = (0, 0),
            ITERATIONS            = 1,
            DRYRUN                = 0,
            #INCLUDE_CHANGED_CELLS = 1,
            SIMULATION_CYKLES     = 0,
            IMAGE_PATH            = 'simulation.result.png'
        )

        if make_settings_file:
            with open("settings.json", "w") as f:
                json.dump(settings, f)

        return cls.from_dict(settings)



    @classmethod
    def from_dict(cls, d: dict):
        return cls(
            P1_PROBABILITY        = float(d[P1_PROBABILITY]),
            P2_PROBABILITY        = float(d[P2_PROBABILITY]),
            DEFECT                = int(d[DEFECT]),
            PERMANENT             = int(d[PERMANENT]),
            FRACTURE_PATH         = int(d[FRACTURE_PATH]),
            DESTRUCTION_PATH      = int(d[DESTRUCTION_PATH]),
            P1                    = int(d[P1]),
            P2                    = int(d[P2]),
            HEAL_CYKLES           = int(d[HEAL_CYKLES]),
            EDGE                  = int(d[EDGE]),
            P1_COLOR              = tuple([int(num) for num in d[P1_COLOR]]),
            P2_COLOR              = tuple([int(num) for num in d[P2_COLOR]]),
            D_COLOR               = tuple([int(num) for num in d[D_COLOR]]),
            PD_COLOR              = tuple([int(num) for num in d[PD_COLOR]]),
            FP_COLOR              = tuple([int(num) for num in d[FP_COLOR]]),
            DP_COLOR              = tuple([int(num) for num in d[DP_COLOR]]),
            BACKGROUND_COLOR      = tuple([int(num) for num in d[BACKGROUND_COLOR]]),
            SUPERCELL_A           = int(d[SUPERCELL_A]),
            SUPERCELL_B           = int(d[SUPERCELL_B]),
            SUPERCELL_C           = int(d[SUPERCELL_C]),
            SUPERCELL_D           = int(d[SUPERCELL_D]),
            SUPERCELL_E           = int(d[SUPERCELL_E]),
            SUPERCELL_F           = int(d[SUPERCELL_F]),
            GRID_SIZE             = tuple([int(num) for num in d[GRID_SIZE]]),
            GRID_SIZE_PIXELS      = tuple([int(num) for num in d[GRID_SIZE_PIXELS]]),
           # Z_AXIS                = int(d[Z_AXIS]),
            FPS                   = int(d[FPS]),
            GRID_MARGINS          = tuple([int(num) for num in d[GRID_MARGINS]]),
            ITERATIONS            = int(d[ITERATIONS]),
            DRYRUN                = int(d[DRYRUN]),
            #INCLUDE_CHANGED_CELLS = int(d[INCLUDE_CHANGED_CELLS]),
            SIMULATION_CYKLES     = int(d[SIMULATION_CYKLES]),
            IMAGE_PATH            = str(d[IMAGE_PATH])
        )

    @classmethod
    def load_from_file(cls, name):
        with open(name, 'r') as f:
            d: dict = json.load(f)

            return cls.from_dict(d)


    def to_dict(self):
        return asdict(self)


    def save_to_file(self, name):
        with open(name, 'w') as f:
            json.dump(asdict(self), f)


    



