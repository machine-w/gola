package gola

import (
	"math"
)

type Matrix3 [9]float64


func (m *Matrix3) SetFromComponents(a, b, c *Vector3) {
	m[0] = a[0]
	m[1] = b[0]
	m[2] = c[0]
	m[3] = a[1]
	m[4] = b[1]
	m[5] = c[1]
	m[6] = a[2]
	m[7] = b[2]
	m[8] = c[2]
}

func (m *Matrix3) Transform(v *Vector3) *Vector3 {
	return &Vector3{
		v[0]*m[0] + v[1]*m[1] + v[2]*m[2],
		v[0]*m[3] + v[1]*m[4] + v[2]*m[5],
		v[0]*m[6] + v[1]*m[7] + v[2]*m[8],
	}
}

func (m *Matrix3) TransformTranspose(v *Vector3) *Vector3 {
	return &Vector3{
		v[0]*m[0] + v[1]*m[3] + v[2]*m[6],
		v[0]*m[1] + v[1]*m[4] + v[2]*m[7],
		v[0]*m[2] + v[1]*m[5] + v[2]*m[8],
	}
}

func (m *Matrix3) TransformMatrix3(b *Matrix3) *Matrix3 {
	newMatrix := &Matrix3{}
	newMatrix[0] = m[0]*b[0] + m[1]*b[3] + m[2] + b[6]
	newMatrix[1] = m[0]*b[1] + m[1]*b[4] + m[2] + b[7]
	newMatrix[2] = m[0]*b[2] + m[1]*b[5] + m[2] + b[8]
	newMatrix[3] = m[3]*b[0] + m[4]*b[3] + m[5] + b[6]
	newMatrix[4] = m[3]*b[1] + m[4]*b[5] + m[5] + b[7]
	newMatrix[5] = m[3]*b[2] + m[4]*b[6] + m[5] + b[8]
	newMatrix[6] = m[6]*b[0] + m[7]*b[3] + m[8] + b[6]
	newMatrix[7] = m[6]*b[1] + m[7]*b[4] + m[8] + b[7]
	newMatrix[8] = m[6]*b[2] + m[7]*b[5] + m[8] + b[8]
	return newMatrix
}

func (m *Matrix3) SetInverse(b *Matrix3) {

	t1 := b[0] * b[4]
	t2 := b[0] * b[5]
	t3 := b[1] * b[3]
	t4 := b[2] * b[3]
	t5 := b[1] * b[6]
	t6 := b[2] * b[6]

	det := t1*b[8] - t2*b[7] - t3*b[8] + t4*b[7] + t5*b[5] - t6*b[4]

	if det == 0 {
		return
	}
	invd := 1 / det

	m[0] = (b[4]*b[8] - b[5]*b[7]) * invd
	m[1] = -(b[1]*b[8] - b[2]*b[7]) * invd
	m[2] = (b[1]*b[5] - b[2]*b[4]) * invd

	m[3] = -(b[3]*b[8] - b[5]*b[6]) * invd
	m[4] = (b[0]*b[8] - t6) * invd
	m[5] = -(t2 - t4) * invd

	m[6] = (b[3]*b[7] - b[4]*b[6]) * invd
	m[7] = -(b[0]*b[7] - t5) * invd
	m[8] = (t1 - t3) * invd

}


func (m *Matrix3) NewInverse() *Matrix3 {
	result := &Matrix3{}
	result.SetInverse(m)
	return result
}

func (m *Matrix3) Invert() {
	m.SetInverse(m)
}

func (m *Matrix3) setInertiaTensorCoeffs(ix, iy, iz, ixy, ixz, iyz float64) {
	m[0] = ix
	m[1] = -ixy
	m[3] = -ixy
	m[2] = -ixz
	m[6] = -ixz
	m[4] = iy
	m[5] = -iyz
	m[7] = -iyz
	m[8] = iz
}

func (m *Matrix3) SetBlockInertiaTensor(halfSizes *Vector3, mass float64) {
	squares := halfSizes.NewHadamardProduct(halfSizes)

	m.setInertiaTensorCoeffs(
		0.3*mass*(squares[1]+squares[2]),
		0.3*mass*(squares[0]+squares[2]),
		0.3*mass*(squares[0]+squares[1]),
		0,
		0,
		0,
	)
}

func (m *Matrix3) SetTranspose(b *Matrix3) {
	m[0] = b[0]
	m[1] = b[3]
	m[2] = b[6]
	m[3] = b[1]
	m[4] = b[4]
	m[5] = b[7]
	m[6] = b[2]
	m[7] = b[5]
	m[8] = b[8]
}

func (m *Matrix3) Transpose() *Matrix3 {
	result := &Matrix3{}
	result.SetTranspose(m)
	return result
}

func (m *Matrix3) SetOrientation(q *Quaternion) {
	m[0] = 1 - (2*q.J*q.J + 2*q.K*q.K)
	m[1] = 2*q.I*q.J + 2*q.K*q.R
	m[2] = 2*q.I*q.K - 2*q.J*q.R
	m[3] = 2*q.I*q.J - 2*q.K*q.R
	m[4] = 1 - (2*q.I*q.I + 2*q.K*q.K)
	m[5] = 2*q.J*q.K + 2*q.I*q.R
	m[6] = 2*q.I*q.K + 2*q.J*q.R
	m[7] = 2*q.J*q.K - 2*q.I*q.R
	m[8] = 1 - (2*q.I*q.I + 2*q.J*q.J)
}

func (m *Matrix3) LinearInterpolate(a, b *Matrix3, prop float64) *Matrix3 {
	result := &Matrix3{}
	for i := uint8(0); i < 9; i++ {
		result[i] = a[i]*(1-prop) + b[i]*prop
	}
	return result
}

type Matrix4 [12]float64

func (m *Matrix4) Clone() *Matrix4 {
	n := &Matrix4{}
	for i := range m {
		n[i] = m[i]
	}
	return n
}

func (m *Matrix4) Equals(other *Matrix4) bool {
	if other == nil {
		return false
	}
	for i := range m {
		if math.Abs(m[i]-other[i]) > RealEpsilon {
			return false
		}
	}
	return true
}

func (m *Matrix4) TransformVector3(b *Vector3) *Vector3 {
	return &Vector3{
		b[0]*m[0] + b[1]*m[1] + b[2]*m[2] + m[3],
		b[0]*m[4] + b[1]*m[5] + b[2]*m[6] + m[7],
		b[0]*m[8] + b[1]*m[9] + b[2]*m[10] + m[11],
	}
}

func (m *Matrix4) TransformMatrix4(b *Matrix4) *Matrix4 {
	newMatrix := &Matrix4{}

	newMatrix[0] = b[0]*m[0] + b[4]*m[1] + b[8]*m[2]
	newMatrix[4] = b[0]*m[4] + b[4]*m[5] + b[8]*m[6]
	newMatrix[8] = b[0]*m[8] + b[4]*m[9] + b[8]*m[10]

	newMatrix[1] = b[1]*m[0] + b[5]*m[1] + b[9]*m[2]
	newMatrix[5] = b[1]*m[4] + b[5]*m[5] + b[9]*m[6]
	newMatrix[9] = b[1]*m[8] + b[5]*m[9] + b[9]*m[10]

	newMatrix[2] = b[2]*m[0] + b[6]*m[1] + b[10]*m[2]
	newMatrix[6] = b[2]*m[4] + b[6]*m[5] + b[10]*m[6]
	newMatrix[10] = b[2]*m[8] + b[6]*m[9] + b[10]*m[10]

	newMatrix[3] = b[3]*m[0] + b[7]*m[1] + b[11]*m[2] + m[3]
	newMatrix[7] = b[3]*m[4] + b[7]*m[5] + b[11]*m[6] + m[7]
	newMatrix[11] = b[3]*m[8] + b[7]*m[9] + b[11]*m[10] + m[11]

	return newMatrix
}

func (m *Matrix4) getDeterminant() float64 {
	return m[8]*m[5]*m[2] + m[4]*m[9]*m[2] + m[8]*m[1]*m[6] - m[0]*m[9]*m[6] - m[4]*m[1]*m[10] + m[0]*m[5]*m[10]
}

func (m *Matrix4) SetInverse(b *Matrix4) {
	det := m.getDeterminant()
	if det == 0 {
		return
	}

	det = 1.0 / det

	m[0] = (-b[9]*b[6] + b[5]*b[10]) * det
	m[4] = (b[8]*b[6] - b[4]*b[10]) * det
	m[8] = (-b[8]*b[5] + b[4]*b[9]) * det

	m[1] = (b[9]*b[2] - b[1]*b[10]) * det
	m[5] = (-b[8]*b[2] + b[0]*b[10]) * det
	m[9] = (b[8]*b[1] - b[0]*b[9]) * det

	m[2] = (-b[5]*b[2] + b[1]*b[6]) * det
	m[6] = (+b[4]*b[2] - b[0]*b[6]) * det
	m[10] = (-b[4]*b[1] + b[0]*b[5]) * det

	m[3] = (+b[9]*b[6]*b[3] - b[5]*b[10]*b[3] - b[9]*b[2]*b[7] + b[1]*b[10]*b[7] + b[5]*b[2]*b[11] - b[1]*b[6]*b[11]) * det
	m[7] = (-b[8]*b[6]*b[3] + b[4]*b[10]*b[3] + b[8]*b[2]*b[7] - b[0]*b[10]*b[7] - b[4]*b[2]*b[11] + b[0]*b[6]*b[11]) * det
	m[11] = (+b[8]*b[6]*b[3] - b[4]*b[9]*b[3] - b[8]*b[1]*b[7] + b[0]*b[9]*b[7] + b[4]*b[1]*b[11] - b[0]*b[5]*b[11]) * det
}

func (m *Matrix4) Inverse() *Matrix4 {
	result := &Matrix4{}
	result.SetInverse(m)
	return result
}

func (m *Matrix4) SetOrientation(q *Quaternion, b *Vector3) {
	m[0] = 1 - (2*q.J*q.J + 2*q.K*q.K)
	m[1] = 2*q.I*q.J + 2*q.K*q.R
	m[2] = 2*q.I*q.K - 2*q.J*q.R
	m[3] = b[0]

	m[4] = 2*q.I*q.J - 2*q.K*q.R
	m[5] = 1 - (2*q.I*q.I + 2*q.K*q.K)
	m[6] = 2*q.J*q.K + 2*q.I*q.R
	m[7] = b[1]

	m[8] = 2*q.I*q.K + 2*q.J*q.R
	m[9] = 2*q.J*q.K - 2*q.I*q.R
	m[10] = 1 - (2*q.I*q.I + 2*q.J*q.J)
	m[11] = b[2]
}


func (m *Matrix4) TransformInverse(b *Vector3) *Vector3 {
	tmp := &Vector3{}
	tmp[0] -= m[3]
	tmp[1] -= m[7]
	tmp[2] -= m[11]

	result := &Vector3{}
	result[0] = tmp[0]*m[0] + tmp[1]*m[4] + tmp[2]*m[8]
	result[1] = tmp[0]*m[1] + tmp[1]*m[5] + tmp[2]*m[9]
	result[2] = tmp[0]*m[2] + tmp[1]*m[6] + tmp[2]*m[10]
	return result
}

func (m *Matrix4) TransformDirection(b *Vector3) *Vector3 {
	result := &Vector3{}
	result[0] = b[0]*m[0] + b[1]*m[1] + b[2]*m[2]
	result[1] = b[0]*m[4] + b[1]*m[5] + b[2]*m[6]
	result[2] = b[0]*m[8] + b[1]*m[9] + b[2]*m[10]
	return result
}

func (m *Matrix4) TransformInverseDirection(b *Vector3) *Vector3 {
	result := &Vector3{}
	result[0] = b[0]*m[0] + b[1]*m[4] + b[2]*m[8]
	result[1] = b[0]*m[1] + b[1]*m[5] + b[2]*m[9]
	result[2] = b[0]*m[2] + b[1]*m[6] + b[2]*m[10]
	return result
}
