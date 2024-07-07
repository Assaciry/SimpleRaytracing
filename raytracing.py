
import pygame as pg
import numpy as np
import json, sys, os

class Sphere():
    def __init__(self, position = (0,0,1), radius = 1, color = (255,0,0)):
        self.position = np.array(position)
        self.radius   = radius
        self.color    = color

    def intersect(self, O, D):
        CO = O - self.position

        a = np.dot(D, D)
        b = 2*np.dot(CO,D)
        c = np.dot(CO,CO) - self.radius**2

        disc = b*b - 4*a*c

        if disc < 0:
            return float("inf"), float("inf")

        t1 = (-b + np.sqrt(disc))/(2*a)
        t2 = (-b - np.sqrt(disc))/(2*a)

        return t1,t2

class Scene():
    def __init__(self, cwidth, cheight, Vw, Vh, f, spheres, BG_COLOR = (255,255,255)):
        self.Cw = cwidth
        self.Ch = cheight
        self.Vw = Vw
        self.Vh = Vh
        self.f  = f

        self.spheres = spheres

        self.BG_COLOR = BG_COLOR

        self.canvas = np.zeros((cheight, cwidth, 3))

    def canvas_to_viewport(self, x,y):
        return x*self.Vw/self.Cw,y*self.Vh/self.Ch, self.f

    def put_pixel(self, x, y, color):
        sx = int(x + self.Cw/2)
        sy = int(y - self.Ch/2)
        self.canvas[sx,sy] = color

    def trace_ray(self, O, D, tmin=1, tmax=float("inf")):
        closest_t = float("inf")
        closest_sphere = None
        for sphere in self.spheres:
            t1,t2 = sphere.intersect(O, D)

            if (t1 >= tmin and t1 <= tmax) and t1 < closest_t:
                closest_t = t1
                closest_sphere = sphere

            if (t2 >= tmin and t2 <= tmax) and t2 < closest_t:
                closest_t = t2
                closest_sphere = sphere

        if closest_sphere == None:
            return self.BG_COLOR

        return closest_sphere.color