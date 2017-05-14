package main

import "math"

type Vector3 struct {
	X, Y, Z float64
}

func (v Vector3) Mul(s float64) Vector3 {
	return Vector3{X: v.X * s, Y: v.Y * s, Z: v.Z * s}
}

func (v Vector3) Div(s float64) Vector3 {
	if s == 0 {
		panic("Div by zero!")
	}
	return Vector3{X: v.X / s, Y: v.Y / s, Z: v.Z / s}
}

func (v Vector3) Add(b Vector3) Vector3 {
	return Vector3{X: v.X + b.X, Y: v.Y + b.Y, Z: v.Z + b.Z}
}

func (v Vector3) Sub(b Vector3) Vector3 {
	return Vector3{X: v.X - b.X, Y: v.Y - b.Y, Z: v.Z - b.Z}
}

func (v Vector3) Inv() Vector3 {
	return Vector3{X: -v.X, Y: -v.Y, Z: -v.Z}
}

func (v Vector3) Dot(b Vector3) float64 {
	return (v.X * b.X) + (v.Y * b.Y) + (v.Z * b.Z)
}

func (v Vector3) MulVec(b Vector3) Vector3 {
	return Vector3{
		X: v.X * b.X,
		Y: v.Y * b.Y,
		Z: v.Z * b.Z,
	}
}

func (v Vector3) Norm() float64 {
	return v.Dot(v)
}

func (v Vector3) Normalize() Vector3 {
	return v.Mul(1.0 / math.Sqrt(v.X*v.X+v.Y*v.Y+v.Z*v.Z))
}

func (v Vector3) CrossProduct(b Vector3) Vector3 {
	return Vector3{
		X: v.Y*b.Z - v.Z*b.Y,
		Y: v.Z*b.X - v.X*b.Z,
		Z: v.X*b.Y - v.Y*b.X,
	}
}
