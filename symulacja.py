import time
import random

from typing import Iterable
from settings import *
from settings import Settings
from unionfind import DisjointSet


import graph

CYKLES_LIMIT_REACHED: int = 100
PATH_FOUND          : int = 101
SIMULATION_RUNNING  : int = 102

class Symulacja2D:
    NS = 0
    WE = 1

    def __init__(self, settings: Settings) -> None:
        '''
        Fields:
        ----------
        __frame: list[list[Cell]]
            list (matrix) of Cells

        __changed: list[tuple[int, int, int]]
            list od cells that have changed since previoues generation

        __width: int
            world width ( number of cells in axis X)

        __height: int
            world height (number of cells in axix Y)

        __datacollector: Datacollector
            for data collection

        __unionfind: UnionFind
            for finding when defects accumulated and destroyed grid


        __found_path: list(tuples[int, int, int])
            contains list of cells that destroyed grid
     
        '''
        self.__frame: list[list[int]] = []
        self.__changed: list[tuple(int, int, int)] = []
        self.settings: Settings = settings

        self._width, self._height  = self.settings.GRID_SIZE

        self._disjoint_set_AC = DisjointSet()
        self._disjoint_set_BD = DisjointSet()
        
        self.edge_cells = set()
        self._seed = self._width
        self._path_starting_point: int = None

        self.make_world()
        self._make_disjoint_sets()

        self._P1_counter: int = self._width * self._height
        self._P2_counter: int = 0
        self._PD_counter: int = 0
        self._D_counter : int = 0

        self.simulation_running: bool = True
        self.cykles_counter: int = 0

        self.connectivity_test_result = None
        self.defecty_d = []
        
    
    def __str__(self):
        frame1: list[list[int]] = self.__frame
        s1 = '\n'.join([''.join(['{:6}'.format(item) for item in row]) for row in frame1])

        return '\n' + s1 + '\n'


    def _encode_id(self, y, x):
        return y* (self._seed + 2) + x


    def _decode_id(self, id_: int):
        x, y = divmod(id_, self._seed + 2)
        return x, y
        

    @staticmethod
    def _vicinity(y:int, x:int) -> tuple[int, int]:
        '''
        Return coordinates of neighbor cell (von Neumann definition)
        
        Parameters:
        ----------
        x: int
            x - axis coordinate of base cell

        y: int
            y - axis coordinate of base cell

        Returns:
        -------
            tuple(int, int) - calculated coordinates of neighbores cells
        '''

        yield (y + 1, x) # Top neighbor
        yield (y - 1, x) # Bottom neighbor
        yield (y, x + 1) # Right neighbor
        yield (y, x - 1) # Left neighbor


  

    def _make_disjoint_sets(self) -> None:
        self.__a = set()
        self.__b = set()
        self.__c = set()
        self.__d = set()
        encode = self._encode_id

        w = self._width
        h = self._height

        unionAC = self._disjoint_set_AC
        unionBD = self._disjoint_set_BD

        unionAC.makeset(self.settings.SUPERCELL_A)
        unionBD.makeset(self.settings.SUPERCELL_B)
        unionAC.makeset(self.settings.SUPERCELL_C)
        unionBD.makeset(self.settings.SUPERCELL_D)

        for i in range(1, w + 1):
            #A
            cell_id = encode(0, i)
            unionAC.makeset(cell_id)
            unionAC.union(cell_id, self.settings.SUPERCELL_A)
            self.__a.add(cell_id)

            #C
            cell_id = encode(h + 1, i) 
            unionAC.makeset(cell_id)
            unionAC.union(cell_id, self.settings.SUPERCELL_C)
            self.__c.add(cell_id)

        for i in range(1, h + 1):
            #B
            cell_id = encode(i, w + 1)
            unionBD.makeset(cell_id)
            unionBD.union(cell_id, self.settings.SUPERCELL_B)
            self.__b.add(cell_id)


            #D
            cell_id = encode(i, 0)
            unionBD.makeset(cell_id)
            unionBD.union(cell_id,self.settings.SUPERCELL_D)
            self.__d.add(cell_id)
            
        self.edge_cells  = self.__a | self.__b | self.__c | self.__d



    @property
    def cykles(self) -> int:
        return self.cykles_counter


    @property
    def PD(self):
        return self._PD_counter

    @property
    def P1(self):
        return self._P1_counter


    @property
    def P2(self):
        return self._P2_counter

    
    @property
    def D(self):
        return self._D_counter

    @property
    def paths_direction(self):
        return self.connectivity_test_result



    def simulation_finished(self) -> bool:
        '''
        Checks if simulation ended. Simulation ends when defects riches oposite edges

        Returns:
        -------
        bool
            True if defects damaged grid otherwise False
        '''
        A = self.settings.SUPERCELL_A
        B = self.settings.SUPERCELL_B
        C = self.settings.SUPERCELL_C
        D = self.settings.SUPERCELL_D

        AC = self._disjoint_set_AC
        BD = self._disjoint_set_BD

        testNS = AC.connected(A, C)
        testWE = BD.connected(B, D)

        if testNS or testWE:
            self.connectivity_test_result = (testNS, testWE)
            if testNS and testWE:
                choice = random.choice([Symulacja2D.NS, Symulacja2D.WE])
                if choice == Symulacja2D.NS:
                    self._path_starting_point = A
                else:
                    self._path_starting_point = B

            elif testWE:
                self._path_starting_point = B

            elif testNS:
                self._path_starting_point = A

            return True

        else:
             return False

       

    def simulation_state(self) -> int:
        # Sprawdź czy znaleziono ścieżke perkolacyjną
        if self.simulation_finished():
            return PATH_FOUND

        # wtakim razie symulacja trwa dalej
        return SIMULATION_RUNNING



    def get_fracture_path(self, mark_connected_nodes: bool = False) -> Iterable[tuple[int, int, int]]:

        start = time.time()
        A: int = self.settings.SUPERCELL_A
        B: int = self.settings.SUPERCELL_B
        C: int = self.settings.SUPERCELL_C
        D: int = self.settings.SUPERCELL_D

        AC: DisjointSet = self._disjoint_set_AC
        BD: DisjointSet = self._disjoint_set_BD

        edge_A: set = self.__a
        edge_B: set = self.__b
        edge_C: set = self.__c
        edge_D: set = self.__d

        supercells = set([A,B,C,D])

        Graph = graph.Graph(weighted=False, digraph=False)
        starting_point: int = self._path_starting_point
        cells = []
        path  = []


        #wez wszystkie weżły które utworzyły rozdarcie w siatce
        if starting_point == A or starting_point == C:
            cells = AC.return_members_set(starting_point)

        elif starting_point == B or starting_point == D:
            cells = BD.return_members_set(starting_point)

        # print("PATHFINDING: Building graph -> start")
        building_graph = time.time()
        for cell in cells:
            Graph.add_node(cell)


        for node in Graph.nodes():
            y, x = self._decode_id(node)
            for b, a in self._vicinity(y, x):

                neighbor = self._encode_id(b, a)
                if neighbor in Graph \
                    and not (neighbor in edge_B and node in edge_D) \
                    and not (neighbor in edge_D and node in edge_B):
                    Graph.add_edge(node, neighbor)

        #połącz superkomórki z komórkami brzegowymi
        if starting_point == A or starting_point == C:
            for node in edge_A:
                if node in Graph:
                    Graph.add_edge(A, node)

            for node in edge_C:
                if node in Graph:
                    Graph.add_edge(C, node)

            path = Graph.shortest_path_bfs(C, A)

        elif starting_point == B or starting_point == D:
            for node in edge_B:
                if node in Graph:
                    Graph.add_edge(B, node)

            for node in edge_D:
                if node in Graph:
                    Graph.add_edge(D, node)

            path = Graph.shortest_path_bfs(D, B)
        
        bg = f'PATHFINDING: Building graph -> end {time.time() - building_graph: 6f}'
        # print(bg)

        returned_elements = []

        #usun superkomorki i komorki brzegowe z listy
        if mark_connected_nodes:
            for elem in cells:
                if elem not in supercells and elem not in self.edge_cells:
                    y, x = self._decode_id(elem)
                    returned_elements.append((y - 1, x - 1, self.settings.DESTRUCTION_PATH))

        for elem in path:
            if elem not in supercells and elem not in self.edge_cells:
                y, x = self._decode_id(elem)
                returned_elements.append((y - 1, x - 1, self.settings.FRACTURE_PATH))

        s = f'PATHFINDING: All {time.time() - start: 6f}'
        # print(s)

        return returned_elements

    

    def make_world(self):
        # print("SIMULATION: BUILDING GRID")
        start = time.time()
        yaxis = self._height + 2
        xaxis = self._width + 2
        frame = self.__frame

        for row in range(yaxis):
            line = []
            for col in range(xaxis):
                if row > 0 and col > 0 and col < self._width + 1 and row < self._height + 1:    
                    line.append(self.settings.P1)

                else:
                    line.append(self.settings.EDGE)

            frame.append(line)

        s = f'\n SUMULATION: GRID BUILDING -> end {time.time() - start}'
        # print(s)
                
               

    def _changed(self, y, x, state) -> None:
        '''
        co

        Parameters:
        ----------
        x: Cell
            Cell object from which we wll extract infomrmation about location
            and category of a cell
        '''

        self.__changed.append(( y - 1, x - 1, state))

    @staticmethod
    def neighbors(y, x):
        yield (y - 1, x)
        yield (y + 1, x)
        yield (y, x - 1)
        yield (y, x + 1)


    def neighbors2(self, y, x, line, next_frame):
        yield (y - 1, x, next_frame[y - 1][x])
        yield (y + 1, x, self.__frame[y + 1][x])
        yield (y, x - 1, line[x - 1])
        yield (y, x + 1, self.__frame[y][x + 1])



    def nerby_defect(self, y, x):
        f = self.__frame
        s = 0
        DEFECT = self.settings.DEFECT

        for b, a in self.neighbors(y, x):
            if f[b][a] < 0:
                s = s + DEFECT
            elif f[b][a] == DEFECT:
                s = s + 0
            else:
                s = s + f[b][a]

        return True if s > DEFECT else False




    def losuj(self, marker_P1, marker_P2):
        frame = self.__frame

        s = self.settings

        for y in range(1, self._height + 1):
            for x in range( 1, self._width + 1):
                cell = self.__frame[y][x]

                if cell == s.P1:
                    if random.random() < s.P1_PROBABILITY:
                        frame[y][x] = marker_P1

                
                elif cell == s.P2:
                    if random.random() < s.P2_PROBABILITY:
                        frame[y][x] = marker_P2

                elif cell > s.DEFECT and cell < s.PERMANENT:
                    frame[y][x] = frame[y][x] - 1
                

  



    def new_cell_state(self, markerp1, markerp2):

        s = self.settings
        next_frame = []

        for y in range(self._height + 2):
            line = []
            for x in range(self._width + 2):
                cell = self.__frame[y][x]
                if x > 0 and y > 0 and x < self._width + 1 and y < self._height + 1:
                    
                    if cell == markerp1 or cell == markerp2:
                        if self.nerby_defect(y, x) == True:
                            cell = s.PERMANENT
                            self._changed(y, x, s.PERMANENT)
                            enc = self._encode_id(y, x)
                            self._disjoint_set_AC.makeset(enc)
                            self._disjoint_set_BD.makeset(enc)

                            if line[x - 1] == s.P1:
                                line[x -1] = s.P2
                                self._changed(y, x - 1, s.P2)
                            
                            if next_frame[y - 1][x] == s.P1:
                                next_frame[y - 1][x] = s.P2
                                self._changed(y - 1, x, s.P2)

                            for b, a, state in self.neighbors2(y, x, line, next_frame):
                                n = state
                                if n == s.PERMANENT or n == s.EDGE:
                                    nid = self._encode_id(b, a)
                                    self._disjoint_set_AC.union(enc, nid)
                                    self._disjoint_set_BD.union(enc, nid)

                            if cell == markerp1:
                                self._P1_counter -= 1
                            else:
                                self._P2_counter -= 1
                            self._PD_counter += 1

                        else:
                            cell = s.DEFECT + s.HEAL_CYKLES
                            if line[x - 1] == s.P1:
                                line[x - 1] = s.P2
                                self._changed(y, x - 1, s.P2)

                            if next_frame[y - 1][x] == s.P1:
                                next_frame[y - 1][x] = s.P2
                                self._changed(y - 1, x, s.P2)
                            
                            # if self.__frame[y][x + 1] == s.P1:
                            #     self.__frame[y][x+1] = s.P2
                            #     self._changed(y, x + 1, s.P2)

                            # if self.__frame[y+1][x] == s.P1:
                            #     self.__frame[y+1][x] = s.P2
                            #     self._changed(y + 1, x, s.P2)
                            
                            self._changed(y, x, s.DEFECT)
                            if cell == markerp1:
                                self._P1_counter -= 1
                            else:
                                self._P2_counter -= 1
                            self._D_counter += 1


                    elif cell == s.P1:
                        if self.nerby_defect(y, x) == True:
                            cell = s.P2
                            self._changed(y, x, cell)
                            self._P1_counter -= 1
                            self._P2_counter += 1

                        else:
                            cell == s.P1

                    elif cell == s.P2:
                        if self.nerby_defect(y, x) == True:
                            cell = s.P2

                        else:
                            cell = s.P1
                            self._changed(y, x, cell)
                            self._P2_counter -= 1
                            self._P1_counter += 1

                    elif cell == s.DEFECT:
                        if line[x - 1] >= s.DEFECT or next_frame[y - 1][x] >= s.DEFECT or \
                            self.__frame[y + 1][x] >= s.DEFECT or self.__frame[y][x + 1] >= s.DEFECT:
                            cell = s.P2
                        else:
                            cell = s.P1
                        self._changed(y, x, cell)
                        self._P1_counter += 1
                        self._D_counter -= 1
                      

                    elif cell > s.DEFECT and cell < s.PERMANENT:
                        if self.nerby_defect(y, x) == True:
                            cell = s.PERMANENT
                            self._changed(y, x, cell)
                            # if line[x - 1] == s.P1:
                            #     line[x - 1] = s.P2
                            #     self._changed(y, x - 1, s.P2)

                            # if next_frame[y-1][x]==s.P1:
                            #     next_frame[y-1][x] = s.P2
                            #     self._changed(y-1, x, s.P2)

                            enc = self._encode_id(y, x)
                            self._disjoint_set_AC.makeset(enc)
                            self._disjoint_set_BD.makeset(enc)

                            for b, a, state in self.neighbors2(y, x, line, next_frame):
                                n = state
                                if n == s.PERMANENT or n == s.EDGE:
                                    nid = self._encode_id(b, a)
                                    self._disjoint_set_AC.union(enc, nid)
                                    self._disjoint_set_BD.union(enc, nid)
                            
                            self._D_counter -= 1
                            self._PD_counter += 1

                        

                line.append(cell)

            else:
                next_frame.append(line)

        return next_frame
                        

                
    def next_step(self):
        self.__changed = []
        marker_p1 = -1
        marker_p2 = -10
        self.losuj(marker_p1, marker_p2)         
        self.__frame = self.new_cell_state(marker_p1, marker_p2)

        self.cykles_counter += 1
            
        return self.__changed


    
    def show_connections(self):
        a = self._disjoint_set_AC.parents
        b = self._disjoint_set_BD.parents

        AC = self._disjoint_set_AC
        BD = self._disjoint_set_BD

        dd = {
            1000020: "B",
            1000010: "A",
            1000030: "C",
            1000040: "D"
        }

        for y in range(self._height+ 2):
            s = f' '
            s2 = f' '
            for x in range(self._width +2):
                elem_id = self._encode_id(y, x)
                if x > 0 and y > 0 and x < self._width + 1 and y < self._height + 1:
                    if elem_id in a:
                        s += ' O '
                    else:
                        s += ' . '
                else:
                    if elem_id in a:
                        s +=' X '
                    else:
                        s += ' | '

            for x in range(self._width +2):
                elem_id = self._encode_id(y, x)
                if x > 0 and y > 0 and x < self._width + 1 and y < self._height + 1:
                    if elem_id in b:
                        s2 += ' O '
                    else:
                        s2 += ' . '
                else:
                    if elem_id in b:
                        s2 +=' X '
                    else:
                        s2 += ' | '
            x = f'{s:1}            {s2:1}'
            # print(x)

        r1 = self._disjoint_set_AC.connected(self.settings.SUPERCELL_A, self.settings.SUPERCELL_C)
        r2 = self._disjoint_set_BD.connected(self.settings.SUPERCELL_B, self.settings.SUPERCELL_D)

        # print('   ',r1,'                                                                    ', r2)


        w = '\n'.join([" ".join(['{:5}'.format(self._encode_id(y,x)) for x in range(self._width + 2)]) for y in range(self._height + 2)])
        
        # # print(w)


        for y in range(self._height + 2):
            s = f' '

            for x in range(self._width + 2):
                elem = self._encode_id(y, x)
                if x > 0 and y > 0 and x < self._width + 1 and y < self._height + 1:
                    if elem in a:
                        p = AC.find(elem)
                        if p in dd:
                            s += f'{dd[p]:3}'
                        else:
                            s += f' {p} '
                        
                    else:
                        s += f' . '

                else:
                    if elem in a:
                        p = AC.find(elem)
                        if p in dd:
                            s += f'{dd[p]:3}'
                        else:
                            s += f' {p} '
                        
                    else:
                        s += f' _ '

            # print(s)

        # print("----------------------------------------------------------")

        for y in range(self._height + 2):
            s = f' '

            for x in range(self._width + 2):
                elem = self._encode_id(y, x)
                if x > 0 and y > 0 and x < self._width + 1 and y < self._height + 1:
                    if elem in b:
                        p = BD.find(elem)
                        if p in dd:
                            s += f'{dd[p]:3}'
                        else:
                            s += f' {p} '
                        
                    else:
                        s += f' . '

                else:
                    if elem in b:
                        p = BD.find(elem)
                        if p in dd:
                            s += f'{dd[p]:3}'
                        else:
                            s += f' {p} '
                        
                    else:
                        s += f' _ '

            # print(s)






        

# class Symulacja3D:
#     """
#     Class that implement 3D naive model of materials destruction
#     It is represent in-memory as a 3D cube of fixed length edge, that 
#     contains cells (like a Rubic's cube). Cube consist of inner cube
#     (1, edge_length + 1), (1, edge_length + 1), (1, edge_length + 1)
#     and outer mantle of 1 layer-deep cells  (0-1, 0-1, 0-1)


#     """



#     def __init__(self, settings: Settings, data_collector = Collector) -> None:
        
#         self.cube = []
#         self.changed = []

#         self.edge_length = settings.Z_AXIS

#         self.path_AC = DisjointSet()
#         self.path_BD = DisjointSet()
#         self.path_EF = DisjointSet()

#         self.graph = graphite.Graph(digraph=False, weighted=False)

#         self.edge_planes = set()
#         self.setings = settings

#         self.P1_counter = self.edge_length ** 3
#         self.P2_counter = 0
#         self.D_counter = 0
#         self.PD_counter = 0
#         self.I_counter = 0
#         self.data_collector = data_collector

#         self.edge_plane_A = set()
#         self.edge_plane_B = set()
#         self.edge_plane_C = set()
#         self.edge_plane_D = set()
#         self.edge_plane_E = set()
#         self.edge_plane_F = set()

#         self.make_world()
#         # print("world created")
#         self.disjoint_sets()
#         # print("disjoint set done")




#     @staticmethod
#     def vicinity( z, y, x):
#         yield(z + 1, y, x)
#         yield(z - 1, y, x)
#         yield(z, y + 1, x)
#         yield(z, y - 1, x)
#         yield(z, y, x + 1)
#         yield(z, y, x - 1)


#     def vicinity2(self, z, y, x, cube, next_cube, plane, line):

#         """
#         Generator that yields coordinates and states of neigbors 


#         We can't depend only on previues cube. We have to have state in next fram from 
#         neighbors that have been already visited. Those negihbors are
#         (z - 1, y, x), (z, y - 1, x) and (z, y, x - 1). With this information we can
#         properly describe next state of current cell

#         Parameters:
#         ===========
#         z, y, x: int
#             cell coordinates in cube
#         cube:
#             previoes cells state
#         next_cube:
#             next cell's states
#         plane:
#             current processed plane in cube with future states of cells
#         line:
#             current processd line in cube with future states of cells

#         Returns:
#         ========
#         Generator: tuple[int, int, int, int]
#             coordinates and states of neigbors. Neighbors:
#             (z - 1, y, x), (z, y - 1, x), (z, y, x -1) will be returned
#             not with current states but with future onees

#         """

#         yield(z - 1, y, x, next_cube[z - 1][y][x]) # Face E site neighbor
#         yield(z + 1, y, x, cube[z + 1][y][x])      # Face F site neighbor
#         yield(z, y - 1, x, plane[y - 1][x])        # Face A site neighbor
#         yield(z, y + 1, x, cube[z][y + 1][x])      # Face C site neighbor
#         yield(z, y, x - 1, line[x - 1])            # Face D site neighbor
#         yield(z, y, x + 1, cube[z][y][x - 1])      # Face B site neighbor

        


#     def save_to_changed(self, z: int, y: int, x: int, state: int) -> None:
#         """
#         Writes changed cells to the changed cells list, and convert 
#         inner cube (1, 1, 1) origin based coodrinates to (0, 0, 0)
#         origin based coordinates. Basically it substract 1 from z, x, y
#         values

#         Parameters:
#         ==========
#         z, x, y: int
#             inner cube z-axis, y-axis, x-axis coordinates

#         state: int
#             cell state
#         """
#         self.changed.append((z - 1, y - 1, x - 1, state))



#     def make_world(self) -> None:
#         """
#         Creates in-memory representation of cube. 
#         All inner-cube cells are P1, and all mantle 
#         cells are EDGE cells
#         """

#         #axes are enlarged by 2 because we have EDGE cells
#         #from two sides on all axes

#         edge_length = self.edge_length + 2
#         cube = self.cube

#         for z in range(edge_length):
#             plane = []
            
#             for y in range(edge_length):
#                 line = []

#                 for x in range(edge_length):
#                     #all inner cells are P1
#                     if  z > 0 and y > 0 and x > 0 \
#                         and z < edge_length - 1 and y < edge_length- 1 \
#                         and x < edge_length - 1:
#                         line.append(self.setings.P1)
#                     #all outer (mantle) cells are EDGE cells
#                     else:
#                         line.append(self.setings.EDGE)

#                 plane.append(line)
#             cube.append(plane)



#     def encode_id(self, z, y, x) -> int:
#         """
#         Encodes cell coordinates to number that
#         is used as a cell ID

#         Parameter:
#         =========
#         z, y, x: int
#             cell coordinates in the cube

#         Returns:
#         ========
#         int:
#             calculated ID
#         """

#         seed = self.edge_length + 2
#         return y * seed + x + z * seed**2



#     def decode_id(self, encoded_id: int) -> tuple[int, int, int]:
#         """
#         Decodes cell ID, to it's (z, y, x)
#         coordinates

#         Parameters:
#         ==========
#         encoded_id: int
#             cell ID

#         Returns:
#         =======
#         tuple[int, int, int]
#             decoded cell coordinates (z, y, x)

#         """

#         seed = self.edge_length + 2
#         z = encoded_id // seed ** 2
#         encoded_id -= z * seed ** 2
#         y, x = divmod(encoded_id, seed)

#         return z, y, x



#     def disjoint_sets(self) -> None:

#         #Locally cache references to 
#         #Disjoint sets
#         pAC = self.path_AC
#         pBD = self.path_BD
#         pEF = self.path_EF

#         #Creates sets for supercells
#         pAC.makeset(self.setings.SUPERCELL_A)
#         pAC.makeset(self.setings.SUPERCELL_C)
#         pBD.makeset(self.setings.SUPERCELL_B)
#         pBD.makeset(self.setings.SUPERCELL_D)
#         pEF.makeset(self.setings.SUPERCELL_E)
#         pEF.makeset(self.setings.SUPERCELL_F)
#         z = self.edge_length
  

#         for i in range(1, z + 1):
#             for j in range(1, z + 1):
#                 #Face E
#                 encoded_id = self.encode_id(0, i, j)
#                 pEF.makeset(encoded_id)
#                 pEF.union(encoded_id, self.setings.SUPERCELL_E)
#                 self.edge_plane_E.add(encoded_id)

#                 #Face F
#                 encoded_id = self.encode_id(z + 1, i, j)
#                 pEF.makeset(encoded_id)
#                 pEF.union(encoded_id, self.setings.SUPERCELL_F)
#                 self.edge_plane_F.add(encoded_id)

#                 #Face D
#                 encoded_id = self.encode_id(j, i, 0)
#                 pBD.makeset(encoded_id)
#                 pBD.union(encoded_id, self.setings.SUPERCELL_D)
#                 self.edge_plane_D.add(encoded_id)

#                 #Face B
#                 encoded_id = self.encode_id(j, i, z + 1)
#                 pBD.makeset(encoded_id)
#                 pBD.union(encoded_id, self.setings.SUPERCELL_B)
#                 self.edge_plane_B.add(encoded_id)

#                 #Face A
#                 encoded_id = self.encode_id(i, 0, j)
#                 pAC.makeset(encoded_id)
#                 pAC.union(encoded_id, self.setings.SUPERCELL_A)
#                 self.edge_plane_A.add(encoded_id)

#                 #Face C
#                 encoded_id = self.encode_id(i, z + 1, j)
#                 pAC.makeset(encoded_id)
#                 pAC.union(encoded_id, self.setings.SUPERCELL_C)
#                 self.edge_plane_C.add(encoded_id)


#     def defected_neighbors(self, cube: list[list[list[int]]], z:int, y:int, x:int):
#         """
#         Check if neighbors are defekted, we have to check
#         if: (z - 1, y, x), (z, y - 1, x) and (z, y, x - 1)
#         are DEFECTS + 1, because if we check only if are defects, then
#         in the next iteration those cells will change to P1 and we will change to PD
#         and we don't want this.

#         Parameters:
#         ==========
#         cube: list[list[list[int]]]
#             in-memory representation of simulated cube
        
#         z, y, x: int
#             coordinates of cell in cube

#         Returns:
#         ========
#         bool:
#             True if there is a defect False otherwise
#         """

#         a = cube[z - 1][y][x] #Already visited in plane before
#         b = cube[z][y - 1][x] #Already visited above
#         c = cube[z][y][x - 1] #Already visited before current (left neighbor)

#         a = self.setings.P1 if a == self.setings.DEFECT + 1 else a
#         b = self.setings.P1 if b == self.setings.DEFECT + 1 else b
#         c = self.setings.P1 if c == self.setings.DEFECT + 1 else c

#         s = a + \
#             cube[z + 1][y][x] + \
#             b + \
#             cube[z][y + 1][x] + \
#             c + \
#             cube[z][y][x + 1]

#         return True if s > self.setings.DEFECT else False


#     def simulation_finished(self) -> bool:
#         """
#         Checks if simulation ended

#         Returns:
#         =======
#         bool:
#             True if simulation end condition is reached
#             False otherwise
#         """

#         A = self.setings.SUPERCELL_A
#         B = self.setings.SUPERCELL_B
#         C = self.setings.SUPERCELL_C
#         D = self.setings.SUPERCELL_D
#         E = self.setings.SUPERCELL_E
#         F = self.setings.SUPERCELL_F

#         AC = self.path_AC
#         BD = self.path_BD
#         EF = self.path_EF

#         if AC.connected(A, C):
#             self.path_starting_point = A
#             return True

#         if BD.connected(B, D):
#             self.path_starting_point = B
#             return True

#         if EF.connected(E, F):
#             self.path_starting_point = E
#             return True

#         return False



#     def get_fracture_path(self, mark_connected_nodes: bool = False):

#         epA = self.edge_plane_A
#         epB = self.edge_plane_B
#         epC = self.edge_plane_C
#         epD = self.edge_plane_D
#         epE = self.edge_plane_E
#         epF = self.edge_plane_F

#         A = self.setings.SUPERCELL_A
#         B = self.setings.SUPERCELL_B
#         C = self.setings.SUPERCELL_C
#         D = self.setings.SUPERCELL_D
#         E = self.setings.SUPERCELL_E
#         F = self.setings.SUPERCELL_F

#         supercells = [A, B, C, D, E, F]

#         AC = self.path_AC
#         BD = self.path_BD
#         EF = self.path_EF

#         s_point = self.path_starting_point

#         graph = graphite.Graph(weighted=False, digraph=False)
#         cells = []
#         path = []

#         match s_point:
#             case A,C:
#                 cells = AC.return_members_set(s_point)
#             case B,D:
#                 cells = BD.return_members_set(s_point)
#             case E,F:
#                 cells = EF.return_members_set(s_point)
            
#         for cell in cells:
#             graph.add_node(cell)

#         for node in graph.nodes():
#             z, y, x = self.decode_id(node)
#             for c, b, a in self.vicinity(z, y, x):
#                 neighbor = self.encode_id(c,b,a)
#                 if neighbor in graph\
#                    and not (neighbor in epB and node in epD)\
#                    and not (neighbor in epD and node in epB):
#                    graph.add_edge(node,neighbor)

#         supernode1 = None
#         supernode2 = None
#         plane1 = None
#         plane2 = None
#         match s_point:
#             case A,C:
#                 supernode1 = A; supernode2 = C
#                 plane1= epA; plane2 = epC
#             case B,D:
#                 supernode2 = D;supernode1 = B
#                 plane1 = epB; plane2=epD
#             case E,F:
#                 supernode1 = E; supernode2 = F
#                 plane1 = epE; plane2 = epF

#         for node in plane1:
#             if node in graph:
#                 graph.add_edge(supernode1, node)

#         for node in plane2:
#             if node in graph:
#                 graph.add_edge(supernode2, node)

#         path = graph.shortest_path_bfs(supernode1, supernode2)

#         returned_elements = []

#         mantle = epA | epB | epC| epD | epE| epF
#         if mark_connected_nodes:
#             for elem in cells:
#                 if elem not in supercells and elem not in mantle:
#                     y, x = self.decode_id(elem)
#                     returned_elements.append((y - 1, x - 1, self.setings.DESTRUCTION_PATH))

#         for elem in path:
#             if elem not in supercells and elem not in mantle:
#                 y, x = self.decode_id(elem)
#                 returned_elements.append((y - 1, x - 1, self.setings.FRACTURE_PATH))


#         return returned_elements



#     def next_step(self) -> list[tuple[int, int, int]]:
#         # print(self.I_counter, self.setings.DEFECT, self.setings.HEAL_CYKLES)
#         zaxis = self.setings.Z_AXIS + 2
#         outer = self.setings.Z_AXIS + 1
#         changed = self.save_to_changed
#         cube = self.cube
#         next_cube: list = []

#         setts = self.setings

#         self.changed = []

#         for z in range(zaxis):
#             plane = []
#             for y in range(zaxis):
#                 line = []
#                 for x in range(zaxis):
#                     cell = cube[z][y][x]

#                     if x > 0 and y > 0 and z > 0 and \
#                        x < outer and y < outer and z < outer:

#                         if cell == setts.P1:
#                             if not self.defected_neighbors(cube, z, y, x):
#                                 if random.random() < self.setings.P1_PROBABILITY:
#                                     cell = setts.DEFECT + setts.HEAL_CYKLES
#                                     changed(z, y, x,setts.DEFECT)
#                                     self.P1_counter -= 1
#                                     self.D_counter += 1

#                                 else:
#                                     cell = setts.P1
#                             else:
#                                 # print("here")
#                                 cell = setts.P2
#                                 changed(z, y, x, cell)
#                                 self.P2_counter += 1
#                                 self.P1_counter -= 1


#                         elif cell == setts.P2:
#                             # print("here")
#                             if not self.defected_neighbors(cube, z, y, x):
#                                 cell = setts.P1
#                                 changed(z, y, x, cell)
#                                 self.P1_counter += 1
#                                 self.P2_counter -= 1
#                             else:
#                                 if random.random() < setts.P2_PROBABILITY:
#                                     cell = setts.PERMANENT
#                                     changed(z, y, x, cell)

#                                     cell_id = self.encode_id(z, y, x)
#                                     self.path_AC.makeset(cell_id)
#                                     self.path_BD.makeset(cell_id)
#                                     self.path_EF.makeset(cell_id)

#                                     for c, b, a, state in self.vicinity2(z, y, x, cube, next_cube, plane, line):
#                                         neighbor = state
#                                         if neighbor == setts.PERMANENT or neighbor == setts.EDGE:
#                                             self.path_AC.union(cell_id, self.encode_id(c, b, a))
#                                             self.path_BD.union(cell_id, self.encode_id(c, b, a))
#                                             self.path_EF.union(cell_id, self.encode_id(c, b, a))

#                                     self.P2_counter -= 1
#                                     self.PD_counter += 1

#                                 else:
#                                     cell = setts.P2


#                         elif cell > setts.DEFECT and cell < setts.PERMANENT:
#                             if not self.defected_neighbors(cube, z, y, x) \
#                                 and next_cube[z - 1][y][x] < setts.PERMANENT \
#                                 and plane[y - 1][x] < setts.PERMANENT \
#                                 and line[x - 1] < setts.PERMANENT:

#                                 cell -= 1
#                                 if cell == setts.DEFECT:
#                                     cell = setts.P1
#                                     changed(z, y, x, cell)
#                                     self.P1_counter += 1
#                                     self.D_counter -= 1

#                             else:
#                                 cell = setts.PERMANENT
#                                 cell_id = self.encode_id(z, y, x)
#                                 self.path_AC.makeset(cell_id)
#                                 self.path_BD.makeset(cell_id)
#                                 self.path_EF.makeset(cell_id)

#                                 for c, b, a, state in self.vicinity2:
#                                     neighbor = state
#                                     if neighbor == setts.PERMANENT or neighbor == setts.EDGE:
#                                         self.path_AC.union(cell_id, self.encode_id(c, b, a))
#                                         self.path_BD.union(cell_id, self.encode_id(c, b, a))
#                                         self.path_EF.union(cell_id, self.encode_id(c, b, a))

#                                 changed(z, y, a, cell)
#                                 self.D_counter -= 1
#                                 self.PD_counter += 1

#                     line.append(cell)

#                 else:
#                     plane.append(line)

#             else:
#                 next_cube.append(plane)

#         else:
#             cube = next_cube
#             self.I_counter += 1


#         if self.setings.INCLUDE_CHANGED_CELLS:
#             self.data_collector.add_record(self.P1_counter, self.P2_counter, self.D_counter,self.PD_counter, self.changed)
#         else:
#             self.data_collector.add_record(self.P1_counter, self.P2_counter, self.D_counter,self.PD_counter)
        

#         # print(self.P1_counter, self.P2_counter, self.D_counter, self.PD_counter)
#         return self.changed


        

            





            





                    
            
        
