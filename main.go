package main

import (
	"fmt"
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

type Plane struct {
	shape  Shape
	normal Vector3
}

type Triangle struct {
	normals [3]Vector3
	shape   Shape
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
	if s == 0 {
		panic("Div by zero!")
	}
	return Vector3{x: v.x / s, y: v.y / s, z: v.z / s}
}

func (v Vector3) Add(b *Vector3) Vector3 {
	return Vector3{x: v.x + b.x, y: v.y + b.y, z: v.z + b.z}
}

func (v Vector3) Sub(b *Vector3) Vector3 {
	return Vector3{x: v.x - b.x, y: v.y - b.y, z: v.z - b.z}
}

func (v Vector3) Inv() Vector3 {
	return Vector3{x: -v.x, y: -v.y, z: -v.z}
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

type LightProperties struct {
	color     Vector3
	intensity float64
}

type PointLight struct {
	properties LightProperties
	position   Vector3
}

type Illuminator interface {
	Illuminate(point Vector3) (float64, Vector3, Vector3)
}

func (l PointLight) Illuminate(point Vector3) (float64, Vector3, Vector3) {
	lightDirection := point.Sub(&l.position)
	r2 := lightDirection.Norm()
	distance := math.Sqrt(r2)
	lightDirAttunated := lightDirection.Div(distance)
	lightIntensity := l.properties.color.Mul(l.properties.intensity / (4 * math.Pi * r2))
	return distance, lightDirAttunated, lightIntensity
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
	normal := phit.Sub(&s.shape.position).Normalize()
	return true, Intersection{point: phit, distance: thit, normal: normal}
}

func (p Plane) Intersect(r *Ray) (bool, Intersection) {
	den := r.direction.Dot(&p.normal)
	if math.Abs(den) < 0.0000001 {
		return false, Intersection{}
	}
	rayToPlane := p.shape.position.Sub(&r.origin)
	t := rayToPlane.Dot(&p.normal) / den
	if t < 0 {
		return false, Intersection{}
	}
	phit := r.origin.Add(&r.direction).Mul(t)
	return true, Intersection{point: phit, distance: t, normal: p.normal}
}

func (t Triangle) Intersect(r *Ray) (bool, Intersection) {
	return true, Intersection{}
}

func CastRay(ray *Ray, objects []Intersectable, lights []Illuminator) Vector3 {
	intersection := Trace(ray, objects, 1000)
	if (intersection == Intersection{}) {
		return Vector3{}
	}
	surfaceColor := Vector3{}
	for _, light := range lights {
		distanceToLight, lightDirection, lightIntensity := light.Illuminate(intersection.point)
		invLightDir := lightDirection.Inv().Normalize()
		lightDirWithBias := invLightDir.Mul(0.001)
		hitPointWithBias := intersection.point.Add(&lightDirWithBias)
		shadowRay := Ray{origin: hitPointWithBias, direction: invLightDir, mint: 0, maxt: 1000}
		lightIntersection := TraceShadowRays(&shadowRay, objects, distanceToLight)
		if (lightIntersection != Intersection{}) {
			continue
		}
		lightContribution := lightIntensity.Mul(math.Max(float64(0), intersection.normal.Dot(&invLightDir)))
		surfaceColor = surfaceColor.Add(&lightContribution)
	}
	clamped := Vector3{x: math.Min(255, surfaceColor.x*255), y: math.Min(255, surfaceColor.y*255), z: math.Min(255, surfaceColor.z*255)}
	return clamped
}

func Trace(ray *Ray, objects []Intersectable, tnear float64) Intersection {
	var foremostObject = Intersection{}
	for _, object := range objects {
		hit, intersection := object.Intersect(ray)
		if hit {
			if intersection.distance < tnear {
				tnear = intersection.distance
				foremostObject = intersection
			}
		}
	}
	return foremostObject
}

func TraceShadowRays(ray *Ray, objects []Intersectable, distanceToLight float64) Intersection {
	for _, object := range objects {
		hit, intersection := object.Intersect(ray)
		if hit && intersection.distance < distanceToLight {
			return intersection
		}
	}
	return Intersection{}
}

func main() {
	fmt.Println("lets trace some butts!")
	// A butt
	s1 := Sphere{radius: 1.5, shape: Shape{position: Vector3{x: 1, y: 0, z: 10}}}
	s2 := Sphere{radius: 1.5, shape: Shape{position: Vector3{x: -1, y: 0, z: 10}}}
	p1 := Plane{normal: Vector3{x: 0, y: 1, z: 0}.Normalize(), shape: Shape{position: Vector3{x: 0, y: -1.5, z: 0}}}

	l1 := PointLight{
		properties: LightProperties{color: Vector3{x: 255, y: 255, z: 255}, intensity: float64(10)},
		position:   Vector3{x: 6, y: 15, z: 7}}

	objects := []Intersectable{s1, s2, p1}
	lights := []Illuminator{l1}

	var width int = 800
	var height int = 600
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
			pixel := CastRay(&cameraRay, objects, lights)
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
