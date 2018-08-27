package gola

import (
	"fmt"
	"math"
)

var RealEpsilon float64

func init() {
	RealEpsilon = 0.00001
}

type Vector3 [3]float64

func (v *Vector3) String() string {
	return fmt.Sprintf("[%0.5f, %0.5f, %0.5f]", v[0], v[1], v[2])
}

var (
	UnitX = Vector3{1, 0, 0}
	UnitY = Vector3{0, 1, 0}
	UnitZ = Vector3{0, 0, 1}
)

func NewVector3(x, y, z float64) *Vector3 {
	e := &Vector3{}
	e[0] = x
	e[1] = y
	e[2] = z
	return e
}

func Zero() *Vector3 {
	return &Vector3{0, 0, 0}
}

func Z() *Vector3 {
	return &Vector3{0, 0, 1}
}

func Y() *Vector3 {
	return &Vector3{0, 1, 0}
}

func X() *Vector3 {
	return &Vector3{1, 0, 0}
}

func (v *Vector3) Clone() *Vector3 {
	return &Vector3{
		v[0],
		v[1],
		v[2],
	}
}

func (v *Vector3) Set(x, y, z float64) {
	v[0] = x
	v[1] = y
	v[2] = z
}

func (v *Vector3) Copy(b *Vector3) {
	v[0] = b[0]
	v[1] = b[1]
	v[2] = b[2]
}

func (v *Vector3) Clear() *Vector3 {
	v[0] = 0
	v[1] = 0
	v[2] = 0
	return v
}

func (v *Vector3) Add(b *Vector3) *Vector3 {
	v[0] += b[0]
	v[1] += b[1]
	v[2] += b[2]
	return v
}

func (v *Vector3) NewAdd(b *Vector3) *Vector3 {
	return &Vector3{
		v[0] + b[0],
		v[1] + b[1],
		v[2] + b[2],
	}
}

func (v *Vector3) Sub(b *Vector3) *Vector3 {
	v[0] -= b[0]
	v[1] -= b[1]
	v[2] -= b[2]
	return v
}

func (v *Vector3) NewSub(b *Vector3) *Vector3 {
	return &Vector3{
		v[0] - b[0],
		v[1] - b[1],
		v[2] - b[2],
	}
}

func (v *Vector3) AddScaledVector(b *Vector3, t float64) *Vector3 {
	if math.IsNaN(t) {
		panic("scale value passed to Vector3.AddScaledVector() is NaN")
	}
	v[0] += b[0] * t
	v[1] += b[1] * t
	v[2] += b[2] * t
	return v
}

func (v *Vector3) Inverse() *Vector3 {
	v[0] = -v[0]
	v[1] = -v[1]
	v[2] = -v[2]
	return v
}

func (v *Vector3) NewInverse() *Vector3 {
	return &Vector3{
		-v[0],
		-v[1],
		-v[2],
	}
}

func (v *Vector3) Length() float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

func (v *Vector3) SquareLength() float64 {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
}

func (v *Vector3) Normalize() *Vector3 {
	length := v.Length()
	if length > 0 {
		v.Scale(1 / length)
	}
	return v
}

func (v *Vector3) Scale(t float64) *Vector3 {
	v[0] *= t
	v[1] *= t
	v[2] *= t
	return v
}

func (v *Vector3) NewScale(t float64) *Vector3 {
	return &Vector3{
		v[0] * t,
		v[1] * t,
		v[2] * t,
	}
}

func (v *Vector3) Dot(b *Vector3) float64 {
	return v[0]*b[0] + v[1]*b[1] + v[2]*b[2]
}

// NewCross aka VectorProduct "%"
func (v *Vector3) NewCross(b *Vector3) *Vector3 {
	return &Vector3{
		v[1]*b[2] - v[2]*b[1],
		v[2]*b[0] - v[0]*b[2],
		v[0]*b[1] - v[1]*b[0],
	}

}


func (v *Vector3) NewVectorProduct(b *Vector3) *Vector3 {
	return v.NewCross(b)
}

func (v *Vector3) ScalarProduct(b *Vector3) float64 {
	return v[0]*b[0] + v[1]*b[1] + v[2]*b[2]
}

func (v *Vector3) HadamardProduct(b *Vector3) *Vector3 {
	v[0] *= b[0]
	v[1] *= b[1]
	v[2] *= b[2]
	return v
}

func (v *Vector3) NewHadamardProduct(b *Vector3) *Vector3 {
	return &Vector3{
		v[0] * b[0],
		v[1] * b[1],
		v[2] * b[2],
	}
}

func (v *Vector3) Equals(b *Vector3) bool {
	diff := math.Abs(v[0] - b[0])
	if diff > RealEpsilon {
		return false
	}
	diff = math.Abs(v[1] - b[1])
	if diff > RealEpsilon {
		return false
	}
	diff = math.Abs(v[2] - b[2])
	return diff < RealEpsilon
}

func (v *Vector3) Rotate(q *Quaternion) *Vector3 {
	num12 := q.I + q.I
	num2 := q.J + q.J
	num := q.K + q.K
	num11 := q.R * num12
	num10 := q.R * num2
	num9 := q.R * num
	num8 := q.I * num12
	num7 := q.I * num2
	num6 := q.I * num
	num5 := q.J * num2
	num4 := q.J * num
	num3 := q.K * num
	num15 := ((v[0] * ((1.0 - num5) - num3)) + (v[1] * (num7 - num9))) + (v[2] * (num6 + num10))
	num14 := ((v[0] * (num7 + num9)) + (v[1] * ((1.0 - num8) - num3))) + (v[2] * (num4 - num11))
	num13 := ((v[0] * (num6 - num10)) + (v[1] * (num4 + num11))) + (v[2] * ((1.0 - num8) - num5))

	v[0] = num15
	v[1] = num14
	v[2] = num13
	return v
}

func (v *Vector3) NewRotate(q *Quaternion) *Vector3 {
	return v.Clone().Rotate(q)
}

