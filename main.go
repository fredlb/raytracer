package main

import (
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"time"
)

type Ray struct {
	origin    Vector3
	direction Vector3
	mint      float64
	maxt      float64
}

// could also hold BBox?
type Shape struct {
	position Vector3
	color    Vector3
	emissive Vector3
	material string
}

type Sphere struct {
	radius float64
	shape  Shape
}

type Intersectable interface {
	Intersect(*Ray) (bool, Intersection)
}

type Intersection struct {
	point    Vector3
	distance float64
	normal   Vector3
	shape    Shape
}

type SphereLight struct {
	sphere   Sphere
	emission Vector3
}

func (s Sphere) Intersect(r *Ray) (bool, Intersection) {
	to := r.origin.Sub(s.shape.position)
	b := to.Dot(r.direction)
	c := to.Dot(to) - s.radius*s.radius
	d := b*b - c
	if d > 0 {
		d = math.Sqrt(d)
		t1 := -b - d
		if t1 > 1e-9 {
			phit := r.origin.Add(r.direction.Mul(t1))
			normal := phit.Sub(s.shape.position)
			normal = normal.Normalize()
			return true, Intersection{point: phit, distance: t1, normal: normal, shape: s.shape}
		}
		t2 := -b + d
		if t2 > 1e-9 {
			phit := r.origin.Add(r.direction.Mul(t2))
			normal := phit.Sub(s.shape.position)
			normal = normal.Normalize()
			return true, Intersection{point: phit, distance: t2, normal: normal, shape: s.shape}
		}
	}
	return false, Intersection{}
}

func Radiance(ray *Ray, objects []Intersectable, lights []SphereLight, depth int, emittance float64, rnd *rand.Rand) Vector3 {
	intersection := Trace(ray, objects, 10000000)
	if depth > 10 {
		return Vector3{X: 0, Y: 0, Z: 0}
	}
	if (intersection == Intersection{}) {
		return Vector3{X: 0, Y: 0, Z: 0}
	}

	orientedSurfaceNormal := Vector3{}
	if intersection.normal.Dot(ray.direction) < 0.0 {
		orientedSurfaceNormal = intersection.normal
	} else {
		orientedSurfaceNormal = intersection.normal.Inv()
	}

	brdfMod := intersection.shape.color
	var p float64
	if brdfMod.X > brdfMod.Y && brdfMod.X > brdfMod.Z {
		p = brdfMod.X
	} else if brdfMod.Y > brdfMod.Z {
		p = brdfMod.Y
	} else {
		p = brdfMod.Z
	}

	if depth+1 > 5 {
		if rnd.Float64() < p*0.9 {
			brdfMod = brdfMod.Mul(0.9 / p)
		} else {
			//Russian roulette
			// TODO: This should take emittance in to account
			return Vector3{X: 0, Y: 0, Z: 0}
		}
	}

	if intersection.shape.material == "DIFF" {
		dr1 := 2 * math.Pi * rand.Float64()
		dr2 := rand.Float64()
		dr2s := math.Sqrt(dr2)
		w := orientedSurfaceNormal
		var u Vector3
		if math.Abs(w.X) > float64(0.1) {
			u = Vector3{X: 0, Y: 1, Z: 0}
		} else {
			u = Vector3{X: 1, Y: 0, Z: 0}
		}
		u = u.CrossProduct(w).Normalize()
		v := w.CrossProduct(u)

		d2 := u.Mul(math.Cos(dr1) * dr2s).Add(v.Mul(math.Sin(dr1) * dr2s)).Add(w.Mul(math.Sqrt(1.0 - dr2)))

		indirectRay := Ray{
			origin:    intersection.point,
			direction: d2.Normalize(),
			mint:      0.0001,
			maxt:      1000,
		}

		em := Vector3{X: 0, Y: 0, Z: 0}
		for _, light := range lights {
			sw := light.sphere.shape.position.Sub(intersection.point)
			su := Vector3{}
			if math.Abs(sw.X) > float64(0.1) {
				su = Vector3{X: 0, Y: 1, Z: 0}
			} else {
				su = Vector3{X: 1, Y: 0, Z: 0}
			}
			su = su.CrossProduct(sw)
			su = su.Normalize()
			sv := sw.CrossProduct(su)

			temp := intersection.point.Sub(light.sphere.shape.position)
			lolDot := temp.Norm()
			cosAMax := math.Sqrt(float64(1) - (light.sphere.radius*light.sphere.radius)/lolDot)
			eps1 := rand.Float64()
			eps2 := rand.Float64()
			cosA := float64(1.0) - eps1 + eps1*cosAMax
			sinA := math.Sqrt(float64(1.0) - cosA*cosA)
			phi := 2 * math.Pi * eps2

			lsu := su.Mul(math.Cos(phi) * sinA)
			lsv := sv.Mul(math.Sin(phi) * sinA)
			lsw := sw.Mul(cosA)
			l := lsu.Add(lsv)
			l = l.Add(lsw)
			l = l.Normalize()

			distance := math.Sqrt(lolDot)
			sr := Ray{origin: intersection.point, direction: l, mint: 0.0001, maxt: 100000}
			shadowRay := Trace(&sr, objects, distance)
			lightRay := Trace(&sr, []Intersectable{light.sphere}, distance)
			if (shadowRay == Intersection{} && lightRay != Intersection{}) {
				omega := 2.0 * math.Pi * (1.0 - cosAMax)
				lightAngle := orientedSurfaceNormal.Dot(l)
				lightContribution := brdfMod.MulVec(light.emission.Mul(lightAngle)).Mul(omega).Div(math.Pi)
				em = em.Add(lightContribution)
			}
		}
		return em.Add(brdfMod.MulVec(Radiance(&indirectRay, objects, lights, depth+1, 0.0, rnd)))
	} else if intersection.shape.material == "SPEC" {
		newRay := Ray{
			origin:    intersection.point,
			direction: ray.direction.Sub(intersection.normal.Mul(2).Mul(intersection.normal.Dot(ray.direction))),
			mint:      0,
			maxt:      1000,
		}
		return brdfMod.MulVec(Radiance(&newRay, objects, lights, depth+1, 0.0, rnd))
	}
	return Vector3{X: 0, Y: 0, Z: 0}
}

func Trace(ray *Ray, objects []Intersectable, tnear float64) Intersection {
	var foremostObject = Intersection{}
	inf := tnear
	for _, object := range objects {
		hit, intersection := object.Intersect(ray)
		if hit {
			if intersection.distance < tnear {
				tnear = intersection.distance
				foremostObject = intersection
			}
		}
	}
	if tnear < inf {
		return foremostObject
	}
	return Intersection{}
}

func Clamp(val float64) float64 {
	if val < 0.0 {
		return 0.0
	} else if val > 1.0 {
		return 1.0
	} else {
		return val
	}
}

func ToInteger(val float64) uint8 {
	return uint8(Clamp(val)*255 + 0.5)
}

func main() {
	fmt.Println("lets trace some butts!")

	// A butt
	//p1 := Plane{normal: Vector3{x: 0, y: 1, z: 0}, shape: Shape{position: Vector3{x: 0, y: -10, z: 0}, color: Vector3{x: 0.75, y: 0.75, z: 0.75}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s2 := Sphere{radius: 5.0, shape: Shape{position: Vector3{X: 5.0, Y: -5.2, Z: 50}, color: Vector3{X: 0.85, Y: 0.85, Z: 0.85}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "SPEC"}}
	s9 := Sphere{radius: 5.0, shape: Shape{position: Vector3{X: -7.0, Y: -5.2, Z: 40}, color: Vector3{X: 0.85, Y: 0.85, Z: 0.85}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "DIFF"}}
	//s1 := Sphere{radius: 5.0, shape: Shape{position: Vector3{X: 5.0, Y: -5.0, Z: 40}, color: Vector3{X: 0.75, Y: 0.75, Z: 0.75}, emissive: Vector3{X: 0, Y: 0, Z: 0}}}
	s3 := Sphere{radius: 100000, shape: Shape{position: Vector3{X: 0, Y: -100010, Z: 10}, color: Vector3{X: 0.75, Y: 0.75, Z: 0.75}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "DIFF"}}
	s4 := Sphere{radius: 100000, shape: Shape{position: Vector3{X: -100015, Y: 0, Z: 10}, color: Vector3{X: 0.75, Y: 0.25, Z: 0.25}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "DIFF"}}
	s6 := Sphere{radius: 100000, shape: Shape{position: Vector3{X: 100015, Y: 0, Z: 10}, color: Vector3{X: 0.25, Y: 0.75, Z: 0.25}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "DIFF"}}
	s5 := Sphere{radius: 100000, shape: Shape{position: Vector3{X: 0, Y: 0, Z: 100070}, color: Vector3{X: 0.75, Y: 0.75, Z: 0.75}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "DIFF"}}
	s8 := Sphere{radius: 100000, shape: Shape{position: Vector3{X: 0, Y: 0, Z: -100050}, color: Vector3{X: 0.75, Y: 0.75, Z: 0.75}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "DIFF"}}
	s7 := Sphere{radius: 100000, shape: Shape{position: Vector3{X: 0, Y: 100055, Z: 10}, color: Vector3{X: 0.75, Y: 0.75, Z: 0.75}, emissive: Vector3{X: 0, Y: 0, Z: 0}, material: "DIFF"}}

	sl1 := SphereLight{
		sphere:   Sphere{radius: 1.2, shape: Shape{position: Vector3{X: 0, Y: 30, Z: 40}, emissive: Vector3{X: 0, Y: 0, Z: 0}}},
		emission: Vector3{X: 550, Y: 550, Z: 550},
	}
	objects := []Intersectable{s2, s4, s5, s6, s7, s3, s8, s9}
	sphereLights := []SphereLight{sl1}

	var width int = 512
	var height int = 384

	spp := 32
	ncpu := 4
	runtime.GOMAXPROCS(ncpu)
	ch := make(chan int, height)
	m := image.NewRGBA(image.Rect(0, 0, width, height))
	black := color.RGBA{0, 0, 0, 255}
	draw.Draw(m, m.Bounds(), &image.Uniform{black}, image.ZP, draw.Src)
	fmt.Println(ncpu)

	camera := LookAt(Vector3{0, 0, 0}, Vector3{0, 0, 100}, Vector3{0, 1, 0}, 35)
	start := time.Now()
	for cpu := 0; cpu < ncpu; cpu++ {
		go func(cpu int) {
			rnd := rand.New(rand.NewSource(time.Now().UnixNano()))
			for i := cpu; i < height; i += ncpu {
				for j := 0; j < width; j++ {
					pixel := Vector3{X: 0, Y: 0, Z: 0}
					for s := 0; s < spp; s++ {
						r1 := rnd.Float64()
						r2 := rnd.Float64()
						ray := camera.CastRay(j, i, width, height, r1, r2)
						pixel = pixel.Add(Radiance(&ray, objects, sphereLights, 0, 1, rnd).Div(float64(spp)))
					}
					pixel = Vector3{X: Clamp(pixel.X), Y: Clamp(pixel.Y), Z: Clamp(pixel.Z)}
					pixelColor := color.RGBA{ToInteger(pixel.X), ToInteger(pixel.Y), ToInteger(pixel.Z), 255}
					m.Set(j, i, pixelColor)
				}
				ch <- 1
			}
		}(cpu)
	}

	for i := 0; i < height; i++ {
		<-ch
	}
	end := time.Now()
	fmt.Println(end.Sub(start))

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
