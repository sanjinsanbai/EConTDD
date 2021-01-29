OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
rz(0.196349540849362) q[3];
u2(0,3.53429173528852) q[2];
u2(0,3.92699081698724) q[1];
rz(6.28318530717959) q[0];
cx q[0], q[1];
rz(6.28318530717959) q[1];
u3(0.785398163397448,1.5707963267949,4.71238898038469) q[0];
cx q[0], q[1];
rz(6.28318530717959) q[1];
u3(-0.785398163397448,1.5707963267949,4.71238898038469) q[0];
cx q[0], q[1];
cx q[1], q[0];
cx q[0], q[1];
cx q[1], q[2];
rz(6.28318530717959) q[2];
u3(0.392699081698724,1.5707963267949,4.71238898038469) q[1];
cx q[1], q[2];
u2(0.785398163397448,3.14159265358979) q[2];
u2(0.392699081698724,3.14159265358979) q[1];
cx q[1], q[2];
h q[1];
h q[2];
cx q[1], q[2];
h q[2];
h q[1];
cx q[1], q[2];
cx q[1], q[0];
tdg q[0];
cx q[1], q[0];
h q[1];
t q[0];
cx q[3], q[4];
h q[3];
h q[4];
cx q[3], q[4];
h q[4];
h q[3];
cx q[3], q[4];
cx q[1], q[2];
cx q[2], q[1];
cx q[1], q[2];
cx q[3], q[4];
cx q[4], q[3];
cx q[3], q[4];
cx q[3], q[1];
rz(-0.196349540849362) q[1];
cx q[3], q[1];
rz(0.392699081698724) q[3];
rz(0.196349540849362) q[1];
cx q[0], q[1];
h q[0];
h q[1];
cx q[0], q[1];
h q[1];
h q[0];
cx q[0], q[1];
cx q[3], q[1];
rz(-0.392699081698724) q[1];
cx q[3], q[1];
t q[3];
rz(0.392699081698724) q[1];
cx q[2], q[1];
h q[2];
h q[1];
cx q[2], q[1];
h q[1];
h q[2];
cx q[2], q[1];
cx q[4], q[3];
h q[4];
h q[3];
cx q[4], q[3];
h q[3];
h q[4];
cx q[4], q[3];
cx q[3], q[4];
cx q[4], q[3];
cx q[3], q[4];
cx q[3], q[1];
tdg q[1];
cx q[3], q[1];
h q[3];
t q[1];
cx q[0], q[2];
cx q[2], q[0];
cx q[0], q[2];
cx q[0], q[1];
cx q[1], q[0];
cx q[0], q[1];
