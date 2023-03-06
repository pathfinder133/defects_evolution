
from typing import Any, Hashable
from settings import *


#custom type hinting
Tdictkey = str | int 


class DisjointSet:
    """
        Implements union-find data structure

        ...
        Attributes
        ----------
        __parent : dict
            contains representants from all sets
        __rank : dict
            contains ranks of all sets roots

        Properties
        ----------
        rank(x: Hashable)
            exposes rank value of a given's element parent
        trees
            returns number of trees in forest
        parents:
            returns parent dictionary
        
        Methods
        -------
        makeset(x)
            creates single node disjoint set, and makes given
            element a root of this set
        find(x)
            returns element parent node (representant)
        union(x, y)
            joins two sets in to one
        connected(x, y)
            checks if x and y are in the same set
    """


    def __init__(self) -> None:
        self.__parent: dict[Hashable, int] = {}
        self.__rank: dict[Hashable, int] = {}
        


    def __str__(self) -> str:
        return str(self.__parent)



    def __contains__(self, x):
        return x in self.__parent


    
    @property
    def rank(self, x) -> int:
        """
        Exposes rank of given element's root
        """

        return self.__rank[x]

    

    @property
    def trees(self) -> int:
        """
        Returns number of trees in 'forest'
        """

        return len(self.__parent)

    

    @property
    def parents(self):
        """
            returns read-only main dict
        """

        return self.__parent



    def exists(self, key):
        return key in self.__parent



    def makeset(self, x: Tdictkey) -> None:
        """Creates new one element disjoint set
            and makes this element a root

            Parameters
            ---------
            x : Hashable
                element from which disjoint set will be created
                must be hashable since it will be used as a dict key
        """

        if x not in self.__parent:
            self.__parent[x] = x
            self.__rank[x]  = 0



    def find(self, x:Tdictkey) -> Any:
        """
            Finds root element of a given one
            Find is implemeted with a path compresion technique

            Parameters
            ----------
            x : Hashable
                element of which root we will be looking for
                must be hashable since we use it as a dict key

            Returns
            -------
            Any
                element that is a representant (root) of given element
                in disjoint set

            Raises
            ------
            KeyError
                If element is not in __parent
         """


        if x != self.__parent[x]:
            self.__parent[x] = self.find(self.__parent[x])
        
        return self.__parent[x]


    
    def union(self, x, y) -> None:
        """
            Merge two disjoint sets
            
            It uses ranks to estimate which tree should become
            subtree of another

            x : Hashable
                element of tree that will be merged with another tree

            y : Hashable
                element of tre that will be merged with another tree
        """
        try:
            rootx: Tdictkey = self.find(x)
            rooty: Tdictkey = self.find(y)

            if rootx != rooty:
                rx: int = self.__rank[rootx]
                ry: int = self.__rank[rooty]

                if ry > rx:
                    self.__parent[rootx] = y
                elif ry < rx:
                    self.__parent[rooty] = x

                else:
                    self.__parent[rootx] = y
                    self.__rank[rooty] += 1

        except KeyError:
            pass


    
    def return_sets(self) -> dict[Tdictkey, list]:

        ret: dict = {}

        for node in self.__parent:
            rootx = self.find(node)

            if rootx not in ret:
                ret[rootx] = [node]
            else:
                ret[rootx] = (node)
                

        return ret 



    def return_members_list(self, x):

        rootx = self.find(x)
        ret = []

        for node in self.__parent:
            if rootx == self.find(node):
                ret.append(node)

        return ret



    def get_parents(self) -> set:
        
        parents = set()
        for parent in self.__parent.values():
            parents.add(parent)

        return parents

    
    
    def return_members_set(self, x: Hashable) -> set:

        rootx = self.find(x)
        ret = set()

        for node in self.__parent:
            if rootx == self.find(node):
                ret.add(node)

        return ret



    def connected(self, x, y) -> bool:
        """
        Checks if two elements are connected (are in the same set)

        Parameters:
        ==========
        x:
            element to find connection
        y:
            element to find connection
        """

        return self.find(x) == self.find(y)


    



