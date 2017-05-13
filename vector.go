package main

type Vector3 struct {
	X, Y, Z float64
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

func (v Vector3) Add(b Vector3) Vector3 {
	return Vector3{x: v.x + b.x, y: v.y + b.y, z: v.z + b.z}
}

func (v Vector3) Sub(b Vector3) Vector3 {
	return Vector3{x: v.x - b.x, y: v.y - b.y, z: v.z - b.z}
}

func (v Vector3) Inv() Vector3 {
	return Vector3{x: -v.x, y: -v.y, z: -v.z}
}

func (v Vector3) Dot(b Vector3) float64 {
	return (v.x * b.x) + (v.y * b.y) + (v.z * b.z)
}

func (v Vector3) MulVec(b Vector3) Vector3 {
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

func (v Vector3) CrossProduct(b Vector3) Vector3 {
	return Vector3{
		x: v.y*b.z - v.z*b.y,
		y: v.z*b.x - v.x*b.z,
		z: v.x*b.y - v.y*b.x,
	}
}
