OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
x q[0];
h q[1];
t q[0];
t q[4];
t q[1];
cx q[4],q[0];
h q[1];
h q[1];
cx q[1],q[4];
cx q[0],q[1];
tdg q[4];
cx q[0],q[4];
tdg q[0];
tdg q[4];
t q[1];
cx q[1],q[4];
cx q[0],q[1];
cx q[4],q[0];
h q[1];
h q[0];
t q[4];
t q[2];
t q[0];
cx q[2],q[4];
cx q[0],q[2];
cx q[4],q[0];
tdg q[2];
cx q[4],q[2];
tdg q[4];
tdg q[2];
t q[0];
h q[7];
h q[7];
cx q[0],q[2];
cx q[4],q[0];
cx q[2],q[4];
h q[0];
h q[5];
t q[3];
h q[1];
h q[1];
t q[2];
t q[5];
cx q[2],q[3];
cx q[5],q[2];
cx q[3],q[5];
tdg q[2];
cx q[3],q[2];
tdg q[3];
tdg q[2];
t q[5];
cx q[5],q[2];
cx q[3],q[5];
cx q[2],q[3];
h q[5];
h q[4];
t q[5];
t q[3];
t q[4];
cx q[3],q[5];
cx q[4],q[3];
cx q[5],q[4];
tdg q[3];
cx q[5],q[3];
h q[4];
h q[4];
tdg q[5];
tdg q[3];
t q[4];
h q[0];
h q[0];
cx q[4],q[3];
cx q[5],q[4];
h q[5];
h q[5];
cx q[3],q[5];
h q[4];
h q[3];
t q[1];
t q[0];
t q[3];
h q[0];
h q[0];
cx q[0],q[1];
cx q[3],q[0];
cx q[1],q[3];
tdg q[0];
cx q[1],q[0];
tdg q[1];
tdg q[0];
t q[3];
cx q[3],q[0];
cx q[1],q[3];
h q[4];
h q[4];
cx q[0],q[1];
h q[3];
h q[4];
t q[5];
t q[3];
t q[4];
cx q[3],q[5];
cx q[4],q[3];
cx q[5],q[4];
tdg q[3];
cx q[5],q[3];
tdg q[5];
tdg q[3];
t q[4];
cx q[4],q[3];
cx q[5],q[4];
cx q[3],q[5];
h q[4];
h q[3];
t q[1];
t q[0];
t q[3];
cx q[0],q[1];
cx q[3],q[0];
cx q[1],q[3];
tdg q[0];
cx q[1],q[0];
tdg q[1];
tdg q[0];
t q[3];
cx q[3],q[0];
cx q[1],q[3];
cx q[0],q[1];
h q[3];
h q[5];
t q[3];
t q[2];
t q[5];
h q[6];
h q[6];
cx q[2],q[3];
cx q[5],q[2];
cx q[3],q[5];
tdg q[2];
cx q[3],q[2];
tdg q[3];
tdg q[2];
t q[5];
cx q[5],q[2];
cx q[3],q[5];
cx q[2],q[3];
h q[5];
h q[5];
h q[5];
h q[4];
t q[5];
h q[5];
h q[5];
t q[3];
t q[4];
cx q[3],q[5];
cx q[4],q[3];
cx q[5],q[4];
tdg q[3];
cx q[5],q[3];
tdg q[5];
h q[2];
h q[2];
tdg q[3];
t q[4];
cx q[4],q[3];
cx q[5],q[4];
cx q[3],q[5];
h q[4];
h q[3];
t q[1];
t q[0];
t q[3];
cx q[0],q[1];
cx q[3],q[0];
cx q[1],q[3];
tdg q[0];
cx q[1],q[0];
tdg q[1];
tdg q[0];
t q[3];
cx q[3],q[0];
cx q[1],q[3];
cx q[0],q[1];
h q[3];
h q[4];
t q[5];
t q[3];
t q[4];
h q[3];
h q[3];
cx q[3],q[5];
cx q[4],q[3];
cx q[5],q[4];
tdg q[3];
h q[7];
h q[7];
cx q[5],q[3];
tdg q[5];
h q[2];
h q[2];
tdg q[3];
t q[4];
cx q[4],q[3];
cx q[5],q[4];
cx q[3],q[5];
h q[4];
h q[3];
t q[1];
t q[0];
t q[3];
h q[3];
h q[3];
cx q[0],q[1];
cx q[3],q[0];
cx q[1],q[3];
tdg q[0];
cx q[1],q[0];
tdg q[1];
tdg q[0];
t q[3];
cx q[3],q[0];
cx q[1],q[3];
cx q[0],q[1];
h q[3];
h q[7];
h q[7];
h q[1];
t q[2];
t q[5];
t q[1];
h q[7];
h q[7];
cx q[5],q[2];
cx q[1],q[5];
cx q[2],q[1];
h q[6];
h q[6];
tdg q[5];
cx q[2],q[5];
h q[2];
h q[2];
tdg q[2];
tdg q[5];
t q[1];
cx q[1],q[5];
cx q[2],q[1];
cx q[5],q[2];
h q[1];
h q[5];
t q[0];
t q[3];
t q[5];
h q[3];
h q[3];
cx q[3],q[0];
cx q[5],q[3];
cx q[0],q[5];
tdg q[3];
h q[3];
h q[3];
cx q[0],q[3];
tdg q[0];
tdg q[3];
t q[5];
cx q[5],q[3];
cx q[0],q[5];
cx q[3],q[0];
h q[6];
h q[6];
h q[5];
h q[1];
t q[2];
t q[5];
t q[1];
cx q[5],q[2];
h q[0];
h q[0];
cx q[1],q[5];
cx q[2],q[1];
tdg q[5];
cx q[2],q[5];
tdg q[2];
tdg q[5];
t q[1];
cx q[1],q[5];
cx q[2],q[1];
cx q[5],q[2];
h q[1];
h q[5];
t q[0];
t q[3];
t q[5];
cx q[3],q[0];
cx q[5],q[3];
cx q[0],q[5];
tdg q[3];
cx q[0],q[3];
tdg q[0];
tdg q[3];
t q[5];
cx q[5],q[3];
cx q[0],q[5];
h q[0];
h q[0];
cx q[3],q[0];
h q[5];
h q[4];
t q[1];
h q[1];
h q[1];
t q[0];
t q[4];
h q[0];
h q[0];
cx q[0],q[1];
cx q[4],q[0];
cx q[1],q[4];
tdg q[0];
cx q[1],q[0];
tdg q[1];
tdg q[0];
t q[4];
cx q[4],q[0];
cx q[1],q[4];
cx q[0],q[1];
h q[4];
