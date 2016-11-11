package main

import (
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"log"
	"math"
	"os"
)

type Vector3 struct {
	x, y, z float64
}

type Ray struct {
	origin    Vector3
	direction Vector3
	mint      float64
	maxt      float64
}

// could also hold BBox?
type Shape struct {
	position Vector3
}

type Sphere struct {
	radius float64
	shape  Shape
}

type Triangle struct {
	normals [3]Vector3
	shape    Shape
}

type Intersectable interface {
	Intersect(*Ray) (bool, Intersection)
}

type Intersection struct {
	point    Vector3
	distance float64
	normal   Vector3
}

func (v Vector3) Mul(s float64) Vector3 {
	return Vector3{x: v.x * s, y: v.y * s, z: v.z * s}
}

func (v Vector3) Div(s float64) Vector3 {
	return Vector3{x: v.x / s, y: v.y / s, z: v.z / s}
}

func (v Vector3) Add(b *Vector3) Vector3 {
	return Vector3{x: v.x + b.x, y: v.y + b.y, z: v.z + b.z}
}

func (v Vector3) Sub(b *Vector3) Vector3 {
	return Vector3{x: v.x - b.x, y: v.y - b.y, z: v.z - b.z}
}

func (v Vector3) Dot(b *Vector3) float64 {
	return (v.x * b.x) + (v.y * b.y) + (v.z * b.z)
}

func (v Vector3) Norm() float64 {
	return v.Dot(&v)
}

func (v Vector3) Normalize() Vector3 {
	if v == (Vector3{0, 0, 0}) {
		return v
	}
	return v.Mul(1 / v.Norm())
}

func QuadricSolveReal(a float64, b float64, c float64) (bool, float64, float64) {
	var discriminant float64 = b*b - float64(4)*a*c
	if discriminant < 0 {
		return false, 0, 0
	}
	var rootDiscriminant = math.Sqrt(discriminant)
	var q float64
	if b < 0 {
		q = float64(-0.5) * (b - rootDiscriminant)
	} else {
		q = float64(-0.5) * (b + rootDiscriminant)
	}

	t1 := q / a
	t2 := c / q
	if t1 > t2 {
		t1, t2 = t2, t1
	}
	return true, t1, t2
}

func (s Sphere) Intersect(r *Ray) (bool, Intersection) {
    rayOrigin := r.origin.Sub(&s.shape.position)
    ray := Ray{origin: rayOrigin, direction: r.direction, mint: r.mint, maxt: r.maxt}
    var A = ray.direction.Dot(&ray.direction)
    var B = 2 * ray.direction.Dot(&ray.origin)
    var C = ray.origin.Dot(&ray.origin) - s.radius*s.radius
    solution, t1, t2 := QuadricSolveReal(A, B, C)
    if !solution {
        return false, Intersection{}
    }
    if t1 > ray.maxt || t2 < ray.mint {
        return false, Intersection{}
    }
    thit := t1
    if t1 < ray.mint {
        thit = t2
        if thit > ray.maxt {
            return false, Intersection{}
        }
    }
    phit := r.origin.Add(&r.direction).Mul(thit)
    normal := phit.Sub(&s.shape.position)
    return true, Intersection{point: phit, distance: thit, normal: normal}
}

func (t Triangle) Intersect(r *Ray) (bool, Intersection) {
	return true, Intersection{}
}

func Trace(ray *Ray, objects []Intersectable) Vector3 {
	var foremostObject = Intersection{}
	var tnearest float64 = 1000
	for _, object := range objects {
		hit, intersection := object.Intersect(ray)
		if hit {
			if intersection.distance < tnearest {
				tnearest = intersection.distance
				foremostObject = intersection
			}
		}
	}
	if (foremostObject == Intersection{}) {
		return Vector3{}
	}
	lightDirection := Vector3{x: 0.5, y: 1, z: -1}.Normalize()
	surfaceColor := Vector3{x: 255, y: 255, z: 255}.Mul(math.Max(float64(0), foremostObject.normal.Dot(&lightDirection)))
	return surfaceColor
}

func main() {
	s1 := Sphere{radius: 1, shape: Shape{position: Vector3{x: 0, y: 0, z: 10}}}
	s2 := Sphere{radius: 1, shape: Shape{position: Vector3{x: 2, y: 2, z: 15}}}
	s3 := Sphere{radius: 1, shape: Shape{position: Vector3{x: -7, y: 1.5, z: 30}}}

	objects := []Intersectable{s1, s2, s3}

	var width int = 640
	var height int = 500
	var invWidth float64 = 1 / float64(width)
	var invHeight float64 = 1 / float64(height)
	var fov float64 = 30
	var aspectratio float64 = float64(width) / float64(height)
	var angle float64 = math.Tan(math.Pi * 0.5 * fov / 180)

	m := image.NewRGBA(image.Rect(0, 0, width, height))
	black := color.RGBA{0, 0, 0, 255}
	draw.Draw(m, m.Bounds(), &image.Uniform{black}, image.ZP, draw.Src)
	for i := 0; i < height; i++ {
		for j := 0; j < width; j++ {
			var xx float64 = (2*((float64(j)+float64(0.5))*invWidth) - 1) * angle * aspectratio
			var yy float64 = (1 - 2*((float64(i)+float64(0.5))*invHeight)) * angle
			cameraRay := Ray{
				origin:    Vector3{x: 0, y: 0, z: 0},
				direction: Vector3{xx, yy, 1}.Normalize(),
				mint:      0, maxt: 1000}
			pixel := Trace(&cameraRay, objects)
			if (pixel != Vector3{}) {
				pixelColor := color.RGBA{uint8(pixel.x), uint8(pixel.y), uint8(pixel.z), 255}
				m.Set(j, i, pixelColor)
			}
		}
	}

	fd, err := os.Create("./output.png")
	if err != nil {
		log.Fatal(err)
	}

	err = png.Encode(fd, m)
	if err != nil {
		log.Fatal(err)
	}

	err = fd.Close()
	if err != nil {
		log.Fatal(err)
	}
}
