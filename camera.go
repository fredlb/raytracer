package main

import (
	"math"
)

type Camera struct {
	p, u, v, w Vector3
	m          float64
}

func LookAt(eye, center, up Vector3, fovy float64) Camera {
	c := Camera{}
	c.p = eye
	c.w = center.Sub(eye).Normalize()
	c.u = up.CrossProduct(c.w).Normalize()
	c.v = c.w.CrossProduct(c.u).Normalize()
	c.m = 1 / math.Tan(fovy*math.Pi/360)
	return c
}

func (c *Camera) CastRay(x, y, w, h int, u, v float64) Ray {
	aspect := float64(w) / float64(h)
	px := ((float64(x)+u-0.5)/(float64(w)-1))*2 - 1
	py := ((float64(y)+v-0.5)/(float64(h)-1))*2 - 1
	d := Vector3{}
	d = d.Add(c.u.Mul(-px * aspect))
	d = d.Add(c.v.Mul(-py))
	d = d.Add(c.w.Mul(c.m))
	d = d.Normalize()
	p := c.p
	return Ray{p, d, 0.0001, 10e5}
}
