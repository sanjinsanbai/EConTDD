from Simulation_with_TDD import *
from TDD.TDD_show import TDD_show
from qiskit.tools.visualization import circuit_drawer
from TDD.TDD import get_node_set, Node, get_int_key
import sys
import os
from cmath import isclose


def eq_of_weight(w1,w2):
    if isclose(w1[0],w2[0]) and isclose(w1[1],w2[1]):
        return 1
    return 0

def eq_of_nodes(self, other):
    """decide whether the two nodes have the same information"""
    if self.key == other.key and self.idx == other.idx and 1 == eq_of_weight(self.out_weight, other.out_weight):
        return True
    else:
        return False

def get_node_list(node,node_list=list()):
    """Only been used when counting the node number of a TDD"""
    node_list.append(node)
    for k in range(2):
        if node.successor[k]:
            node_list = get_node_list(node.successor[k],node_list)
    return node_list


def pattern_check(node):
    """to determine which pattern the node v belongs to"""
    '''pattern:
    0: not 1/2/3/4
    1: w0(v)=0 and w0(high(v))=0
    2: w1(v)=0 and w1(low(v))=0
    3: w0(v)=w1(v)=1 and w1(low(v))=0 and w0(high(v))=0 and high(low(v))=low(high(v))
    4: w0(v)=w1(v)=1 and w1(low(v))=0 and w0(high(v))=0 and high(low(v))≠low(high(v))
    '''
    if node.key == 1 or node.key[0] == 'y':
        return 0
    if get_int_key(node.out_weight[0]) == (0, 0):
        if get_int_key(node.successor[1].out_weight[0]) != (0, 0):
            return 0
        else:
            return 1
    if get_int_key(node.out_weight[1]) == (0, 0):
        if get_int_key(node.successor[0].out_weight[1]) != (0, 0):
            return 0
        else:
            return 2
    if node.successor[0].key == node.successor[1].key:
        if get_int_key(node.successor[0].out_weight[1]) != (0, 0) or \
                get_int_key(node.successor[1].out_weight[0]) != (0, 0):
            return 0
        else:
            if get_int_key(node.out_weight[0]) == get_int_key(node.out_weight[1]) and \
                    node.successor[0].successor[0] == node.successor[1].successor[1]:
                return 3
            else:
                return 4
    return 0



def RR5(tdd, cur):
    """to return a tdd that removes redundant nodes and preserves circuit equivalence"""
    '''preorder traversal on the graph'''
   
    if cur.key == 1: 
        return cur
    pattern = pattern_check(cur)
    if pattern == 1:
        cur.key = cur.successor[1].key
        cur.successor[1] = cur.successor[1].successor[1]
        cur.out_weight[1] = 1
    if pattern == 2:
        cur.key = cur.successor[0].key
        cur.successor[0] = cur.successor[0].successor[0]
        cur.out_weight[0] = 1
    if pattern == 3:
        cur =  cur.successor[0].successor[0]
        cur = RR5(tdd,cur)
        return cur
    if pattern == 4:
        cur.key = cur.successor[0].key
        cur.successor[0] = cur.successor[0].successor[0]
        cur.successor[1] = cur.successor[1].successor[1]

    for k in range(2):
        cur.successor[k] = RR5(tdd,cur.successor[k])
    return cur



def equivalence_check(tdd1,tdd2):
    """decide whether the two tdds are isomorphic"""
    node_list1 = list()
    node_list2 = list()
    node_list1 = get_node_list(tdd1.node,node_list1)
    node_list2 = get_node_list(tdd2.node,node_list2)
    len1 = len(node_list1)
    len2 = len(node_list2)
    #print('Node number of tdd1 is ', len1)
    #print('Node number of tdd2 is ', len2)
    if len1 != len2:
        #print("eq-check: reason: len1 != len2")
        return False
    for i in range(len1):
        if not eq_of_nodes(node_list1[i],node_list2[i]):
            #print("eq-check: err at node :",node_list1[i].key,node_list2[i].key)
            #print(node_list1[i].idx,node_list2[i].idx)
            #print(node_list1[i].out_weight,node_list2[i].out_weight)
            return False
    return True



def RR5_and_eqCheck(partionScheme, origin_path, modified_path, file_name):
    '''
        param:
            partion scheme = { 0, 1, 2 }
            origin/modified path ： the path of origin curcuits and modified curcuits
            file_name : 
    '''
    qbits_num = 0
    gate_num = 0

    file_list = []
    file_list.append(file_name)
    file_list.append("out_"+file_name)
    tdd = [None]*2
    flag = 0
    
   
    for file_name in file_list:

        #print("-------------------------------------------------------------------")
        if flag == 0:
            #print("origin",end=" ")
            path = origin_path
        else:
            #print("modified",end=" ")
            path = modified_path

        DGStartTime = time.time()
        cir,res = CreateDGfromQASMfile(file_name, path, flag_single=True)
        #print("DG time: %.6f"%(time.time() - DGStartTime))
        dag_cir=res[0]

        qbits_num = get_real_qubit_num(dag_cir)
        
        gate_num = get_gates_number(dag_cir)
        

        simulationStartTime = time.time()
        tmpTDD,max_node_num,block_num=Simulation_with_cir_partion(dag_cir, qbits_num,partionScheme,qbits_num//2,qbits_num//2+1)
        #print("simulation with cir partion time: %.6f"%(time.time() - simulationStartTime))
        simulation_time = time.time() - simulationStartTime
        tdd[flag] = tmpTDD
        #run_time=time.time()-t_start
        #print('Run time:',run_time)
        #print('Max node num: ',max_node_num)
        final_node_num = tmpTDD.node_number()

        rr5StartTime = time.time()                                           #RR5
        tdd[flag].node = RR5(tdd[flag],tdd[flag].node)                 
        rr5_time = time.time() - rr5StartTime

        # after_rr5_node_num = tdd[flag].node_number()

        rr5_time = -1
        after_rr5_node_num = -1

        #tag show cir and tdd
        # TDD_show(tdd[flag],file_name+"_out")
        # print(cir.draw())
        print(qbits_num,"\t",gate_num,"\t",simulation_time, \
            "\t",max_node_num,"\t",final_node_num,"\t",rr5_time,"\t", \
            after_rr5_node_num,end = "\t")
        
        flag+=1
    #print('=================')
  
    #eqCheck
    #print("-------------------------------------------------------------------")
    t_start = time.time()
    iseq = equivalence_check(tdd[0],tdd[1])
    t_end = time.time()
    print(iseq,"\t",(t_end - t_start))
    #print("equivalence: ",iseq)
    #print("check time: %.6f"%(t_end - t_start))


def batch_check(partionScheme):
    '''
    check the equivalence of curcuits in batch

    param：
        partion scheme = { 0, 1, 2 }
    '''

    origin_path = "D:\\Projects\\Algorithm\\tdd\\benchmark\\origin\\"
    target_path = "D:\\Projects\\Algorithm\\tdd\\benchmark\\modify\\"

    for root, subdir, files in os.walk(origin_path):
        for origin_file in files:

            #print("\n!!===================================================================")
            #print("#circuits name: ",origin_file)
            print(origin_file,end="\t")
            t_start = time.time()
            RR5_and_eqCheck(partionScheme,origin_path,target_path,origin_file)
            t_end = time.time()
            
            #print("TDD time: %.6f"%(t_end - t_start))

def simulation(path, file_name, partionScheme):

        cir,res = CreateDGfromQASMfile(file_name, path, flag_single=True)
        #print("DG time: %.6f"%(time.time() - DGStartTime))
        dag_cir=res[0]

        qbits_num = get_real_qubit_num(dag_cir)
        gate_num = get_gates_number(dag_cir)
        
        simulationStartTime = time.time()
        tmpTDD,max_node_num,block_num=Simulation_with_cir_partion(dag_cir, qbits_num,partionScheme,qbits_num//2,qbits_num//2+1)
        #print("simulation with cir partion time: %.6f"%(time.time() - simulationStartTime))
        simulation_time = time.time() - simulationStartTime
        
        #print('Run time:',run_time)
        #print('Max node num: ',max_node_num)
        final_node_num = tmpTDD.node_number()

        TDD_show(tmpTDD,file_name+"_befor_rr5")

        rr5StartTime = time.time()
        tmpTDD.node = RR5(tmpTDD,tmpTDD.node)                  #RR5
        rr5_time = time.time() - rr5StartTime

        # after_rr5_node_num = tmpTDD.node_number()
        TDD_show(tmpTDD,file_name+"_after_rr5")
        print(cir.draw())
