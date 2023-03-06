 
import math
import pickle
import time
from   multiprocessing import Pipe, Queue
import pygame
from   settings import *
import symulacja

from symulacja import Symulacja2D as Sim2D
from   ui import Grid



data = []
runs_times = []

def run_with_visualisation(settings: Settings):

    pygame.init()
    screen = pygame.display.set_mode(settings.GRID_SIZE_PIXELS)
    clock  = pygame.time.Clock()
    grid   = Grid(settings)
    world  = symulacja.Symulacja2D(settings)

    screen.blit(grid.get_surface(), settings.GRID_MARGINS)
    pygame.display.update()

    cykles_limited : bool = True if settings.SIMULATION_CYKLES > 0 else False
    cykles = settings.SIMULATION_CYKLES

    running = True
    pause = False
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    pause = not pause

        if not pause:
            gen = world.next_step()
            grid.update(gen, screen)

            state = world.simulation_state()
            if state == symulacja.PATH_FOUND:
                grid.update(world.get_fracture_path(mark_connected_nodes=True), screen)
                running = False
            
            elif cykles_limited == True:
                cykles -= 1
                if cykles == 0:
                    running = False
                    
            pygame.display.update()
            clock.tick(settings.FPS)

    timee = f'{time.time()}'

    pygame.image.save(screen, 'simulation_result_' + timee+'.png')
    pygame.quit()

    return (world.P1, world.P2, world.D, world.PD, world.cykles, None, None)



def run_without_visualisation(settings: Settings):
    running = True
    world = symulacja.Symulacja2D(settings)

    cykles_limited : bool = True if settings.SIMULATION_CYKLES > 0 else False
    cykles = settings.SIMULATION_CYKLES

    while running:
        world.next_step()

        state = world.simulation_state()
        if state == symulacja.PATH_FOUND:
            running = False

        elif cykles_limited == True:
            cykles -= 1
            if cykles == 0:
                running = False

    
    paths = world.paths_direction
    if not paths:
        paths = (None, None)

    return (world.P1, world.P2, world.D, world.PD, world.cykles, paths[Sim2D.NS], paths[Sim2D.WE])              
    
    
def  run_simulation(settings: Settings, pipe: Pipe, queue:Queue):

    s = settings
    time.sleep(0.5) # żeby mieć pewność że wszystko wystartowało

    for run in range(s.ITERATIONS):
        start = time.time()
        run_info = []

        simulation_info = (
            run,
            s.P1_PROBABILITY,
            s.P2_PROBABILITY,
            s.HEAL_CYKLES,
            s.GRID_SIZE[0],
            s.GRID_SIZE[1],
            s.SIMULATION_CYKLES
        )
        run_info.append(simulation_info)

        if s.DRYRUN: # bez wizualizacji
            result = run_without_visualisation(settings)

        elif not s.DRYRUN:# z wizualizają
            result = run_with_visualisation(settings)

        run_info.append(result)

        percent    = math.floor(run*100 / s.ITERATIONS + .5)
        cykle_time = time.time() - start
        fps        = 1 / cykle_time
        eta        = (s.ITERATIONS - run) / fps

        run_info.append((cykle_time, eta, fps))
        data.append(run_info)
        queue.put((run, s.ITERATIONS, percent, fps, eta), block=False)

    else:
        pickle.dump(data, open("pickle_last_run", "wb"))
        pipe.send(data)
        pipe.close()






        


    








