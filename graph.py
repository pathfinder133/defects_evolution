import heapq

from collections import defaultdict, deque
from typing import Any, Hashable, Iterator

from unionfind import DisjointSet, Tdictkey


class Graph:
    """
    Implementacja grafu za pomocą słownika jeżyka Python. Klasa oprócz struktury
    zawiera również kilka najprostszych operacji na grafie jak: dodawanie i odejmowanie
    węzłów i krawędzi, trawersowanie w głąb (DFS) i wszerz (BFS), odnajdywanie najkrótszej
    ścieżki (BFS i Dijikstra) a takżę minimalnego drzewa rozpinającego (Kruskal)
    ....
    
    Methods
    -------
    digraph(): property
        returns true if graph is directed and false if not
    weighted(): property
        returns true if graph is weighted and false if not
    num_of_nodes(): property
        returns nnumber of nodes in graph
    graph_name(): property
        returns string contained graph name
        default: empty
    add_node(node)

    """

    def __init__(self, name: str = '', digraph: bool = False, weighted: bool = False) -> None:
        
        self.__graph: dict[Hashable, Any] = {}

        self.__nodes: int = 0
        self.__edges: int = 0

        self.__name: str = name

        #Czy graf jest kierunkowy
        self.__digraph: bool = digraph
        #Czy graf jest ważony
        self.__weighted: bool = weighted



    def __str__(self) -> str:
        """
        Zwraca reprezentacje grafu w formie łańcucha znaków
        """

        return str(self.__graph)


    def __contains__(self, node):
        return True if node in self.__graph else False



    def __len__(self):
        return len(self.__graph)


    
    def __getitem__(self, item):
        return self.__graph[item]



    def __setitem__(self, key, value):
        self.__graph[key] = value



    @property
    def digraph(self) -> bool:
        return self.__digraph



    @property
    def weighted(self) -> bool:
        return self.__weighted



    @property
    def num_of_nodes(self) -> int:
        return self.__nodes

    

    @property
    def graph_name(self) -> str:
        return self.__name


    
    @property
    def num_of_edges(self) -> int:
        x = self.__edges
        return x

    

    @graph_name.setter
    def graph_name(self, name: str) -> None:
        self.__name = name



    def nodes(self) -> Iterator:
        for node in self.__graph.keys():
            yield node

    

    def add_node(self, node) -> None:
        '''
        Dodaje węzeł do grafu
        '''

        #jeżeli graf jest ważony to inaczej trzeba rerezentować połączenia
        #z sąsiadami. Dla ważonego reprezentacją jest słownik w formie
        #k:v = nazwa_sąsiada: waga_polączenia
        #dla nieważonego to tylko nazwa sąsiada
        match self.weighted:
            case True:
                self.__graph[node] = {}

            case False:
                self.__graph[node] = set()

        
        self.__nodes += 1



    def add_edge(self, nodeA: Hashable, nodeB: Hashable, weight: int = None) -> None:
        """
        Dodaje krawędź pomiędzy węzłąmi. Jeżeli węzły nieistnieją
        #to najpierw są tworzone
        """

        graph: dict = self.__graph

        if nodeA not in graph:
            self.add_node(nodeA)

        if nodeB not in graph:
            self.add_node(nodeB)

        match self.digraph, self.weighted:
            case True, True:
                graph[nodeA][nodeB] = weight
            
            case True, False:
                graph[nodeA].add(nodeB)

            case False, True:
                graph[nodeA][nodeB] = weight
                graph[nodeB][nodeA] = weight

            case False, False:
                graph[nodeA].add(nodeB)
                graph[nodeB].add(nodeA)

        self.__edges += 2 if nodeA == nodeB else 1
    


    def get_adjacency(self, node) -> dict | set:
        """
        Zwraca wszystkich sąsiadów węzła. Dla grafu kierunkowego
        zwraca sąsiadóœ na których wskazuję ten węzeł
        """

        try:
            neighbors = self.__graph[node]

        except KeyError:
            raise KeyError

        else:
            return neighbors


    
    def remove_node(self, removenNode) -> None:
        """
        Usuwa węzeł z grafu. Usuwa róznież krawędzie łączące go z innymi węzłąmi
        """

        graph: dict = self.__graph
        try:
            match self.digraph, self.weighted:
                case True, True:
                    for node in graph:
                        if removenNode in graph[node]:
                            del graph[node][removenNode]
                    
                    del graph[removenNode]

                case False, True:
                    for neighbor in self.get_adjacency(removenNode):
                        del graph[neighbor][removenNode]

                    del graph[removenNode]

                case True, False:
                    for node in graph:
                        if removenNode in graph[node]:
                            graph[node].remove(removenNode)

                    del graph[removenNode]

                case False, False:
                    for neighbor in self.get_adjacency(removenNode):
                        graph[neighbor].remove(removenNode)

                    del graph[removenNode]
                
        except KeyError:
            raise KeyError

        
        else:
            self.__nodes -= 1



    def remove_edge(self, nodeA, nodeB) -> None:
        """
        Usuwa krawędź pomiędzy węzłami
        """

        graph: dict = self.__graph
        try:
            match self.digraph, self.weighted:
                case True, True:
                    del graph[nodeA][nodeB]

                case False, True:
                    del graph[nodeA][nodeB]
                    del graph[nodeB][nodeA]

                case True, False:
                    graph[nodeA].remove(nodeB)
                
                case False, False:
                    graph[nodeA].remove(nodeB)
                    graph[nodeB].remove(nodeA)

        except KeyError:
            raise KeyError



    def transpose(self) -> dict:
        """
        Transpozycja grafu kierunkowego
        """

        if self.digraph:
            new_graph = {}

            for node in self.__graph:
                for neighbor in self.get_adjacency():
                    new_graph[neighbor] = node

            return new_graph

        else:
            raise ValueError



    def bfs(self, root) -> Iterator:
        """
        Przejście grafu wszerz.
        """
        
        visited = set()
        queue = deque()

        visited.add(root)
        queue.append(root)

        while queue:
            node = queue.popleft()
            yield node

            for neighbor in self.get_adjacency(node):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)



    def dfs(self, root) -> Iterator:
        """
        Trawers grafu w głąb
        """

        visited: set = set()
        stack: list = []

        visited.add(root)
        stack.append(root)

        while stack:

            node = stack.pop()
            yield node

            for neighbor in self.get_adjacency(node):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)


    
    def shortest_path_bfs(self, nodeA, nodeB) -> list:
        """
        Wyszukuje najkrótszą ścieżke w grafie (w sensie najmniejszej liczby węzłów)
        pomiędzy dwoma węzłami
        za pomocą algorytmu BFS
        """

        visited = set()
        queue = deque()
        path = {}
        found_path = []

        visited.add(nodeA)
        queue.append(nodeA)
        path[nodeA] = None
        success = False

        while queue:
            node = queue.popleft()
            if node == nodeB:
                success = True
                break

            for neighbor in self.get_adjacency(node):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

                    path[neighbor] = node

        if success:
            found_path.append(nodeB)

            next_node = nodeB
            while path[next_node] is not None:
                found_path.append(path[next_node])
                next_node = path[next_node]

        found_path.reverse()


        del visited
        del path
        del queue
        return found_path
                


    def dijkstra(self, nodeA, nodeB):
        """
        Wyszukiwanie najkrótszej ścieżki (w sensie najmniejszej wagi krawędzi)
        za pomocą algorytmu Dijikstry. Implementacja za pomocą 
        kolejki priorytetowej
        """

        graph: dict = self.__graph

        #visited nodes
        visited: set = set()

        #map where at the end, each node will have
        #as a value, node from where dijikstra 
        #traversed to it
        parentsMap: dict = {}

        #priority queue. Will keep, as it root, node with
        #current lowest cost
        pq: list = []

        #map where will be kept cost to all nodes from the starting node
        nodeCost: dict = defaultdict(lambda: float('inf'))

        #travel cost to the starting node to the starting node
        #is of course 0
        nodeCost[nodeA] = 0

        heapq.heappush(pq, (0, nodeA))

        while pq:

            #we don't need cost from priority queue
            #we only need name of the nearest node
            __, node = heapq.heappop(pq)

            visited.add(node)

            for neighbor, distance in graph[node].items():

                if neighbor in visited:
                    continue

                newcost = nodeCost[node] + distance
                
                if nodeCost[neighbor] > newcost:
                    parentsMap[neighbor] = node
                    nodeCost[neighbor] = newcost
                    heapq.heappush(pq, (newcost, neighbor))

        path = []
        curNode = nodeB
        while curNode is not nodeA:
            path.append(curNode)
            curNode = parentsMap[curNode]

        path.append(nodeA)
        path.reverse()

        
        return parentsMap, nodeCost, path



    def SCM(self) -> dict[Tdictkey, list[Tdictkey]]:
        """
        Wyszukiwanie silnie połączonych skłądowych w grafie
        .....
        
        Returns
        -------
            T: dict[Tdictkey, list[Tdictkey]
                dict with all groups of nodes in graph
                dict[group_representant: list[all nodes belonged to group]]
            """

        T = DisjointSet()

        graph = self.__graph

        for node in graph:
            T.makeset(node)

        for node in graph:
            for neighbor in graph[node]:
                T.union(node, neighbor)

        return T.return_sets()



    def get_random_node(self):
        """
        Zwraca losowo węzeł z grafu
        """

        item = self.__graph.popitem()

        self.__graph[item[0]] = item[1]

        return item



    def MST(self):
        """
        Minimalne drzewo rozpinające za pomocą algorytmu Kruskala
        """

        queue = []
        T = []
        Z = DisjointSet()

        graph = self.__graph

        for node in graph:
            Z.makeset(node)

            for neighbor, value in graph[node].items():
                heapq.heappush(queue, (value, neighbor, node))

        while queue:
            _, u, v = heapq.heappop(queue)
            if Z.find(u) == Z.find(v):
                continue
            else:
                T.append((u,v))
                Z.union(u,v)

        return T


    
    def is_cyclic(self) -> bool:
        """
        Sprawdza czy graf jest cykliczny
        """

        visited = set()
        stack = []

        start = self.get_random_node()[0]
        stack.append(start)
        stack.append(-1)

        visited.add(start)

        while stack:

            w = stack.pop()
            v = stack.pop()
            for z in self.get_adjacency(v):
                if z not in visited:
                    stack.append(z)
                    stack.append(v)
                    visited.add(z)

                elif z != w:
                    return True
                
        
        return False

    


