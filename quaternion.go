package gola

import (
	"math"
)
type Quaternion struct {
	R float64
	I float64
	J float64
	K float64
}

func NewQuaternion(r, i, j, k float64) *Quaternion {
	return &Quaternion{r, i, j, k}
}

func QuaternionToTarget(origin, target *Vector3) *Quaternion {
	dest := target.NewSub(origin).Normalize()

	source := Z()
	dot := source.Dot(dest)
	if math.Abs(dot-(-1.0)) < RealEpsilon {
		return QuaternionFromAxisAngle(Y(), -math.Pi)
	} else if math.Abs(dot-(1.0)) < RealEpsilon {
		return &Quaternion{1, 0, 0, 0}
	}
	rotAngle := math.Acos(dot)
	rotAxis := source.NewCross(dest).Normalize()
	return QuaternionFromAxisAngle(rotAxis, rotAngle)
}

func QuaternionFromAxisAngle(axis *Vector3, angle float64) *Quaternion {
	halfSin := math.Sin(angle / 2)
	halfCos := math.Cos(angle / 2)
	q := &Quaternion{
		halfCos,
		axis[0] * halfSin,
		axis[1] * halfSin,
		axis[2] * halfSin,
	}
	return q
}

func QuaternionFromVectors(a, b *Vector3) *Quaternion {

	m := math.Sqrt(2.0 + 2.0*a.Dot(b))

	w := a.Clone().NewCross(b).Scale(1.0 / m)

	return &Quaternion{
		0.5 * m,
		w[0],
		w[1],
		w[2],
	}
}

func (q *Quaternion) Set(r, i, j, k float64) {
	q.R = r
	q.I = i
	q.J = j
	q.K = k
}

func (q *Quaternion) Clone() *Quaternion {
	return &Quaternion{
		R: q.R,
		I: q.I,
		J: q.J,
		K: q.K,
	}
}

func (q *Quaternion) Equals(z *Quaternion) bool {
	if math.Abs(q.R-z.R) > RealEpsilon {
		return false
	}
	if math.Abs(q.I-z.I) > RealEpsilon {
		return false
	}
	if math.Abs(q.J-z.J) > RealEpsilon {
		return false
	}
	if math.Abs(q.K-z.K) > RealEpsilon {
		return false
	}
	return true
}

func (q *Quaternion) Normalize() {
	d := q.R*q.R + q.I*q.I + q.J*q.J + q.K*q.K

	if d < RealEpsilon {
		q.R = 1
		return
	}
	d = 1.0 / math.Sqrt(d)
	q.R *= d
	q.I *= d
	q.J *= d
	q.K *= d
}


func (q *Quaternion) Conjugate() *Quaternion {
	q.I = -q.I
	q.J = -q.J
	q.K = -q.K
	return q
}

func (q *Quaternion) NewConjugate() *Quaternion {
	return q.Clone().Conjugate()
}

func (q *Quaternion) NewInverse() *Quaternion {
	t := q.SquareLength()
	return &Quaternion{
		q.R / t,
		-q.I / t,
		-q.J / t,
		-q.K / t,
	}
}

func (q *Quaternion) Dot(q2 *Quaternion) float64 {
	return q.R*q2.R + q.I*q2.I + q.J*q2.J + q.K*q2.K
}

func (q *Quaternion) Div(s float64) *Quaternion {
	return &Quaternion{q.R / s, q.I / s, q.J / s, q.K / s}
}

func (q *Quaternion) Length() float64 {
	d := q.R*q.R + q.I*q.I + q.J*q.J + q.K*q.K
	return math.Sqrt(d)
}

func (q *Quaternion) SquareLength() float64 {
	return q.R*q.R + q.I*q.I + q.J*q.J + q.K*q.K
}

func (q *Quaternion) Norm() float64 {
	return q.SquareLength()
}

func (q *Quaternion) Multiply(o *Quaternion) *Quaternion {
	*q = Quaternion{
		-q.I*o.I - q.J*o.J - q.K*o.K + q.R*o.R,
		q.I*o.R + q.J*o.K - q.K*o.J + q.R*o.I,
		-q.I*o.K + q.J*o.R + q.K*o.I + q.R*o.J,
		q.I*o.J - q.J*o.I + q.K*o.R + q.R*o.K,
	}
	return q
}

func (q *Quaternion) NewMultiply(o *Quaternion) *Quaternion {
	return q.Clone().Multiply(o)
}

func (q *Quaternion) AddScaledVector(vector *Vector3, scale float64) {

	vectorQ := &Quaternion{0, vector[0] * scale, vector[1] * scale, vector[2] * scale}
	result := vectorQ.NewMultiply(q)
	result.Div(2)

	q.R += result.R
	q.I += result.I
	q.J += result.J
	q.K += result.K
}

func (q *Quaternion) RotateByVector(vector *Vector3) *Quaternion {
	return q.Multiply(&Quaternion{0, vector[0], vector[1], vector[2]})
}

func (q *Quaternion) NewRotateByVector(vector *Vector3) *Quaternion {
	return q.NewMultiply(&Quaternion{0, vector[0], vector[1], vector[2]})
}

func LocalToWorld(local *Vector3, transform *Matrix4) *Vector3 {
	return transform.TransformVector3(local)
}

func WorldToLocal(world *Vector3, transform *Matrix4) *Vector3 {
	return transform.TransformInverse(world)
}

func LocalToWorldDirn(local *Vector3, transform *Matrix4) *Vector3 {
	return transform.TransformDirection(local)
}

func WorldToLocalDirn(world *Vector3, transform *Matrix4) *Vector3 {
	return transform.TransformInverseDirection(world)
}