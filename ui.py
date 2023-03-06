import math
import pygame
from settings import *


class Grid:
    """
    Creates a surface and paints a grid on it. Also update the grid
    according to the list of changes.


    Attributes:
    ==========
    settings: Settings
        dataclass cthat contains all program parameters and settings
    _grid: list[list[pygame.Rect]]
        memory representation of a grid, for easy tranlation row, col coordinates
        to actual coordinates on screen
    _surface: pygame.Surface
        canvas for grid drawing and updating
    _surface_margins: tuple(int,int)
        on screen margins of grid (px)
    _grid_size_in_pixels: tuple(int,int)
        size of grid in pixel

    Methods:
    =======
    _make_grid(self)
        calculating grid parameters, and in-memory grid representation
    get_surface(self)
        returns canvas
    update(cells, surface)
        updates grid in memory and on screen, according to information
        from cells

    
    Imports:
    =======
    pygame
    """


    def __init__(self, settings: Settings) -> None:
        """
        Parameter:
        =========
        settings: Settings
            dataclass that contains all program settings

        """

        self.settings             :Settings                = settings
        self._grid                :list[list[pygame.Rect]] = []
        self._surface             :pygame.Surface          = pygame.Surface(settings.GRID_SIZE_PIXELS)
        self._grid_size           :tuple[int, int]         = settings.GRID_SIZE
        self._surface_margins     :tuple[int, int]         = settings.GRID_MARGINS
        self._grid_size_in_pixels :tuple[int, int]         = settings.GRID_SIZE_PIXELS

        self._surface.fill(self.settings.BACKGROUND_COLOR)
        self._make_grid()


    def _make_grid(self) -> None:
        """
        Calculates coordinates and sizes of grid's tiles
        """

        grid_width, grid_height = self._grid_size_in_pixels
        x, y = self._grid_size
        tile_x = grid_width // x
        tile_y = grid_height // y

        #To avoid tiles being a float number
        #calculates  tiles size as an integer
        #and cumulative difference between all tiles size
        #and screen width is the margin
        margin_x = math.floor((grid_width - tile_x * x) / 2)
        margin_y = math.floor((grid_height - tile_y * y) / 2)

        for j in range(y):
            line = []
            for i in range(x):
                #we need those rect to later quicly change them 
                #without calculation, Classic space for speed trade
                rect = pygame.Rect(margin_x + i * tile_x,
                                   margin_y + j * tile_y,
                                   tile_x - 1,
                                   tile_y - 1)

                line.append(rect)
                pygame.draw.rect(self._surface, color=self.settings.P1_COLOR, rect=rect)
   
            self._grid.append(line)



    def get_surface(self) -> pygame.Surface:
        return self._surface


    
    def update(self, cells:list[tuple[int,int,int]], surface: pygame.Surface) -> None:
        """
        Updates grid with information from changed tiles list

        Parameters:
        ==========
        cells: list[tuple[int,int,int]]
            list of changed tails, consisted from y, x and cell's state
            y, x  = row and col of tile
            state determine color of a tile

        surface: pygame.Surface
            a surface on which tiles will be drawn 
            TODO: why this is not a constructor parameter
                  available through self?
        """

        for y, x, state in cells:

            if state == self.settings.P1:
                    color = self.settings.P1_COLOR

            elif state == self.settings.P2:
                    color = self.settings.P2_COLOR
                
            elif state == self.settings.DEFECT:
                    color =  self.settings.D_COLOR

            elif state == self.settings.PERMANENT:
                    color = self.settings.PD_COLOR

            elif state ==self.settings.DESTRUCTION_PATH:
                    color = self.settings.DP_COLOR

            elif state == self.settings.FRACTURE_PATH:
                    color = self.settings.FP_COLOR

            rect = self._grid[y][x]
            pygame.draw.rect(self._surface, color, rect = rect)

        surface.blit(self._surface, self._surface_margins)
    