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
	color    Vector3
	emissive Vector3
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
	shape    Shape
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

func (v Vector3) Add2(b Vector3) Vector3 {
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

func (v Vector3) MulVector3(b *Vector3) Vector3 {
	return Vector3{
		x: v.x * b.x,
		y: v.y * b.y,
		z: v.z * b.z,
	}
}

func (v Vector3) MulVector32(b Vector3) Vector3 {
	return Vector3{
		x: v.x * b.x,
		y: v.y * b.y,
		z: v.z * b.z,
	}
}

func (v Vector3) Norm() float64 {
	return v.Dot(&v)
}

func (v Vector3) Normalize() Vector3 {
	return v.Mul(1.0 / math.Sqrt(v.x*v.x+v.y*v.y+v.z*v.z))
}

func (v Vector3) CrossProduct(b *Vector3) Vector3 {
	return Vector3{
		x: v.y*b.z - v.z*b.y,
		y: v.z*b.x - v.x*b.z,
		z: v.x*b.y - v.y*b.x,
	}
}

type SphereLight struct {
	sphere   Sphere
	emission Vector3
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
	return true, Intersection{point: phit, distance: thit, normal: normal, shape: s.shape}
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
	return true, Intersection{point: phit, distance: t, normal: p.normal, shape: p.shape}
}

func (t Triangle) Intersect(r *Ray) (bool, Intersection) {
	return true, Intersection{}
}

func Radiance(ray *Ray, objects []Intersectable, lights []SphereLight, depth int, emittance float64, rnd *rand.Rand) Vector3 {
	intersection := Trace(ray, objects, 1000)
	if depth > 10 {
		return Vector3{x: 0, y: 0, z: 0}
	}
	if (intersection == Intersection{}) {
		return Vector3{x: 0, y: 0, z: 0}
	}

	orientedSurfaceNormal := Vector3{}
	if intersection.normal.Dot(&ray.direction) < 0.0 {
		orientedSurfaceNormal = intersection.normal.Normalize()
	} else {
		orientedSurfaceNormal = intersection.normal.Normalize().Inv()
	}

	brdfMod := intersection.shape.color
	var p float64
	if brdfMod.x > brdfMod.y && brdfMod.x > brdfMod.z {
		p = brdfMod.x
	} else if brdfMod.y > brdfMod.z {
		p = brdfMod.y
	} else {
		p = brdfMod.z
	}

	if depth+1 > 5 {
		if rnd.Float64() < p {
			brdfMod = brdfMod.Mul(1.0 / p)
		} else {
			//Russian roulette
			return Vector3{x: 0, y: 0, z: 0}
		}
	}

	r1 := rnd.Float64()
	r2 := rnd.Float64()
	sample := UniformSampleHemisphere(r1, r2)
	Nt, Nb := CreateCoordinateSystem(orientedSurfaceNormal)
	sampleWorld := Vector3{
		x: sample.x*Nb.x + sample.y*intersection.normal.x + sample.z*Nt.x,
		y: sample.x*Nb.y + sample.y*intersection.normal.y + sample.z*Nt.y,
		z: sample.x*Nb.z + sample.y*intersection.normal.z + sample.z*Nt.z,
	}
	indirectRay := Ray{
		origin:    intersection.point.Add(&sampleWorld).Mul(0.0001),
		direction: sampleWorld.Normalize(),
		mint:      0,
		maxt:      1000,
	}
	normalWBias := orientedSurfaceNormal.Mul(0.0001)
	hitPointWithBias := intersection.point.Add(&normalWBias)

	em := Vector3{x: 0, y: 0, z: 0}
	for _, light := range lights {
		sw := light.sphere.shape.position.Sub(&intersection.point)
		su := Vector3{}
		if math.Abs(sw.x) > float64(0.1) {
			su = Vector3{x: 0, y: 1, z: 0}
		} else {
			su = Vector3{x: 1, y: 0, z: 0}
		}
		su = su.CrossProduct(&sw)
		su = su.Normalize()
		sv := sw.CrossProduct(&su)

		temp := hitPointWithBias.Sub(&light.sphere.shape.position)
		lolDot := temp.Norm()
		cosAMax := math.Sqrt(float64(1) - light.sphere.radius*light.sphere.radius/lolDot)
		eps1 := rnd.Float64()
		eps2 := rnd.Float64()
		cosA := float64(1.0) - eps1 + eps1*cosAMax
		sinA := math.Sqrt(float64(1.0) - cosA*cosA)
		phi := 2 * math.Pi * eps2

		lsu := su.Mul(math.Cos(phi) * sinA)
		lsv := sv.Mul(math.Sin(phi) * sinA)
		lsw := sw.Mul(cosA)
		l := lsu.Add(&lsv)
		l = l.Add2(lsw).Normalize()

		distance := math.Sqrt(lolDot)
		sr := Ray{origin: hitPointWithBias, direction: l, mint: 0, maxt: 1000}
		shadowRay := Trace(&sr, objects, distance)
		lightRay := Trace(&sr, []Intersectable{light.sphere}, 1000)
		if (shadowRay == Intersection{} && lightRay != Intersection{}) {
			omega := 2.0 * math.Pi * (1.0 - cosAMax)
			lightAngle := orientedSurfaceNormal.Dot(&l)
			lightContribution := brdfMod.MulVector32(light.emission.Mul(lightAngle * omega)).Div(math.Pi)
			em = em.Add2(lightContribution)
		}
	}
	return em.Add2(brdfMod.MulVector32(Radiance(&indirectRay, objects, lights, depth+1, 0.0, rnd)))
}

func UniformSampleHemisphere(r1 float64, r2 float64) Vector3 {
	sinTheta := math.Sqrt(1.0 - r1*r2)
	phi := 2.0 * math.Pi * r2
	x := sinTheta * math.Cos(phi)
	z := sinTheta * math.Sin(phi)
	return Vector3{x: x, y: r1, z: z}.Normalize()
}

func CreateCoordinateSystem(normal Vector3) (Vector3, Vector3) {
	Nt := Vector3{}
	if math.Abs(normal.x) > math.Abs(normal.y) {
		Nt = Vector3{x: normal.z, y: 0, z: -normal.x}.Div(math.Sqrt(normal.x*normal.x + normal.z*normal.z))
	} else {
		Nt = Vector3{x: 0, y: -normal.z, z: normal.y}.Div(math.Sqrt(normal.y*normal.y + normal.z*normal.z))
	}

	Nb := normal.CrossProduct(&Nt)
	return Nt, Nb
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
	return uint8(math.Pow(Clamp(val), 1.0/2.2)*255 + 0.5)
}

func main() {
	fmt.Println("lets trace some butts!")

	// A butt
	s2 := Sphere{radius: 4.0, shape: Shape{position: Vector3{x: -6.0, y: -6.0, z: 37}, color: Vector3{x: 0.95, y: 0.95, z: 0.95}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s1 := Sphere{radius: 4.0, shape: Shape{position: Vector3{x: 6.0, y: -6.0, z: 30}, color: Vector3{x: 0.95, y: 0.95, z: 0.95}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s3 := Sphere{radius: 100000, shape: Shape{position: Vector3{x: 0, y: -100010, z: 10}, color: Vector3{x: 0.75, y: 0.75, z: 0.75}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s4 := Sphere{radius: 100000, shape: Shape{position: Vector3{x: -100013, y: 0, z: 10}, color: Vector3{x: 0.95, y: 0.25, z: 0.25}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s6 := Sphere{radius: 100000, shape: Shape{position: Vector3{x: 100013, y: 0, z: 10}, color: Vector3{x: 0.25, y: 0.25, z: 0.95}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s5 := Sphere{radius: 100000, shape: Shape{position: Vector3{x: 0, y: 0, z: 100050}, color: Vector3{x: 0.75, y: 0.75, z: 0.75}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s8 := Sphere{radius: 100000, shape: Shape{position: Vector3{x: 0, y: 0, z: -100045}, color: Vector3{x: 0.0, y: 0.0, z: 0.0}, emissive: Vector3{x: 0, y: 0, z: 0}}}
	s7 := Sphere{radius: 100000, shape: Shape{position: Vector3{x: 0, y: 100010, z: 10}, color: Vector3{x: 0.75, y: 0.75, z: 0.75}, emissive: Vector3{x: 0, y: 0, z: 0}}}

	sl1 := SphereLight{
		sphere:   Sphere{radius: 2.0, shape: Shape{position: Vector3{x: 8, y: 4, z: 33}, emissive: Vector3{x: 0, y: 0, z: 0}}},
		emission: Vector3{x: 20, y: 20, z: 20},
	}
	objects := []Intersectable{s3, s1, s2, s4, s5, s6, s7, s8}
	sphereLights := []SphereLight{sl1}

	var width int = 512
	var height int = 384
	var invWidth float64 = 1 / float64(width)
	var invHeight float64 = 1 / float64(height)
	var fov float64 = 40
	var aspectratio float64 = float64(width) / float64(height)
	var angle float64 = math.Tan(math.Pi * 0.5 * fov / 180)

	spp := 256
	ncpu := 4
	runtime.GOMAXPROCS(ncpu)
	ch := make(chan int, height)
	m := image.NewRGBA(image.Rect(0, 0, width, height))
	black := color.RGBA{0, 0, 0, 255}
	draw.Draw(m, m.Bounds(), &image.Uniform{black}, image.ZP, draw.Src)
	fmt.Println(ncpu)
	start := time.Now()
	for cpu := 0; cpu < ncpu; cpu++ {
		go func(cpu int) {
			rnd := rand.New(rand.NewSource(time.Now().UnixNano()))
			for i := cpu; i < height; i += ncpu {
				for j := 0; j < width; j++ {
					pixel := Vector3{x: 0, y: 0, z: 0}
					for s := 0; s < spp; s++ {
						r1 := rnd.Float64()
						r2 := rnd.Float64()
						var xx float64 = (2*((float64(j)+r1)*invWidth) - 1) * angle * aspectratio
						var yy float64 = (1 - 2*((float64(i)+r2)*invHeight)) * angle
						cameraRay := Ray{
							origin:    Vector3{x: 0, y: 0, z: 0},
							direction: Vector3{xx, yy, 1}.Normalize(),
							mint:      0, maxt: 100,
						}

						pixel = pixel.Add2(Radiance(&cameraRay, objects, sphereLights, 0, 1, rnd).Div(float64(spp)))
					}
					pixel = Vector3{x: Clamp(pixel.x), y: Clamp(pixel.y), z: Clamp(pixel.z)}
					pixelColor := color.RGBA{ToInteger(pixel.x), ToInteger(pixel.y), ToInteger(pixel.z), 255}
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
