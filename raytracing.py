import numpy as np
import json, sys, os

from enum import Enum

class LightType(Enum):
    AMBIENT = 0
    POINT   = 1
    DIRECTIONAL = 2

class Light():
    def __init__(self, intensity, position = None, direction = None):
        self.intensity = intensity

        if position is not None and direction is not None:
            raise ValueError("Both position and direction parameters are given. You can only give either position for Point Light or direction for Directional Light")

        elif position is not None:
            self.type = LightType.POINT
            self.position = np.array(position,dtype=np.float32)
        
        elif direction is not None:
            self.type = LightType.DIRECTIONAL
            self.direction = np.array(direction,dtype=np.float32)
        
        else:
            self.type = LightType.AMBIENT

class Sphere():
    def __init__(self, position = (0,0,1), radius = 1, color = (255,0,0), specular_coeff = 100, reflectance = 0.):
        self.position = np.array(position,dtype=np.float32)
        self.radius   = radius
        self.color    = np.array(color,dtype=np.float32)
        self.specular_coeff = specular_coeff
        self.reflectance = reflectance

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
    def __init__(self, cwidth, cheight, Vw, Vh, f, spheres, lights, BG_COLOR = (255,255,255)):
        self.Cw = cwidth
        self.Ch = cheight
        self.Vw = Vw
        self.Vh = Vh
        self.f  = f

        self.spheres = spheres
        self.lights  = lights
        self.BG_COLOR = np.array(BG_COLOR,dtype=np.float32)
        self.canvas = np.zeros((cheight, cwidth, 3))

    def canvas_to_viewport(self, x,y):
        return x*self.Vw/self.Cw,y*self.Vh/self.Ch, self.f

    def put_pixel(self, x, y, color):
        color = color / 255.
        sx = int(x + self.Cw/2)
        sy = int(y - self.Ch/2)
        self.canvas[sx,sy] = color

    def trace_ray(self, O, D, tmin=1, tmax=float("inf"), recursion_depth=0):
        closest_sphere, closest_t = self.closest_intersection(O,D,tmin,tmax)

        if closest_sphere == None:
            return self.BG_COLOR

        # Add lighting
        P = O + closest_t * D
        N = P - closest_sphere.position
        N = N / np.linalg.norm(N)
        light_intensity = self.compute_lighting(P, N, -D, closest_sphere.specular_coeff)
        local_color = closest_sphere.color * light_intensity

        closest_reflectance = closest_sphere.reflectance
        if recursion_depth <= 0 or closest_reflectance <= 0:
            return local_color
        
        R = self.reflect_ray(-D, N)
        reflected_color = self.trace_ray(P, R, 0.001, float("inf"), recursion_depth-1)

        return closest_reflectance * reflected_color + (1 - closest_reflectance) * local_color

    @staticmethod
    def reflect_ray(R, N):
        return 2*N*np.dot(N,R)-R

    def closest_intersection(self, O, D, tmin=1, tmax=float("inf")):
        """ Calculate the ray-sphere intersection for a ray starting from _O_ along the direction _D_. Returns _(sphere,t)_"""
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

        return closest_sphere, closest_t

    def compute_lighting(self, p, n, v, s):
        i = 0. # Light intensity
        
        for light in self.lights:
            if light.type == LightType.AMBIENT:
                i += light.intensity
            else:
                if light.type == LightType.POINT:
                    L = light.position - p
                    tmax = 1
                elif light.type == LightType.DIRECTIONAL:
                    L = light.direction
                    tmax = float("inf")

                # Shadow check
                shadow_sphere, shadow_t = self.closest_intersection(p, L, 0.001, tmax)
                if shadow_sphere is not None:
                    continue

                # Diffuse light
                n_dot_L = np.dot(n, L)
                if n_dot_L > 0:
                    i += light.intensity * n_dot_L/(np.linalg.norm(n)*np.linalg.norm(L))
                
                # Specular light
                if s != -1: # set s = -1 if the object is matte
                    R = 2*n*n_dot_L-L
                    r_dot_v = np.dot(R,v)
                    if r_dot_v > 0:
                        rn = np.linalg.norm(R)
                        vn = np.linalg.norm(v)
                        i_ = i
                        i += light.intensity * (r_dot_v/(rn*vn))**s 

        return i