from IndexElimination import batch_check, simulation, RR5_and_eqCheck


if __name__ == '__main__':

    partionScheme = 1
#eq checking
#batch
    # print("#eqCheck running...TIME:",end = ' ')
    # print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) )
    # batch_check(partionScheme)

#single
    # t_start = time.time()
    # RR5_and_eqCheck(partionScheme),"D:\\Projects\\Algorithm\\tdd\\benchmark\\","D:\\Projects\\Algorithm\\tdd\\benchmark\\","test.qasm")
    # t_end = time.time()
    # print("TDD time:",t_end - t_start)
  
#simulation
    simulation("./benchmark/","test.qasm",partionScheme)