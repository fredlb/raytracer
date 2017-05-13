package main

import "testing"

func TestAddVector(t *testing.T) {
	vec1 := Vector3{X: 1, Y: 2, Z: 3}
	vec2 := Vector3{X: 1, Y: 2, Z: 3}
	expected := Vector3{X: 2, Y: 4, Z: 6}
	res := vec1.Add(vec2)
	if res != expected {
		t.Error("Expected", expected, "got", res)
	}
}
