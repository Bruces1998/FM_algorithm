##### IMPORTS ######
import random
import datetime
import pickle as pl
import os

#### GLOBAL VARIABLES ####
sizeA = 0
sizeB = 0
sizeA_minimum_cut = 0
sizeB_minimum_cut = 0
cutset = 0
new_cutset = 0
temp_index = 0  
locked_cell_index = 0
minimum_cut = 0
startcut = 0
old_partition = 0

locked_cells = []
gain_bucket = {}
max_gain_index = None  # index for gainA
best_cut_array = []
best_cut = float("inf")
best_copies_cell_map = []

cell_map = {}
net_map = {}

class Cell:
    def __init__(self, id):
        self.id = id
        self.cell_size = 0
        self.F = 0
        self.T = 0
        self.cell_gain = 0
        self.lock_status = 0  # 1=locked, 0=unlocked
        self.cell_type = ""
        self.cell_partition = 0  # A = 0, B = -1
        self.net_list = set()
        self.cell_size = 0

class Net:
    def __init__(self, id):
        self.id = id
        self.cell_list = []
        self.cutstate = 0  # 0 = uncut, -1 = cut
        self.Asize = 0
        self.Bsize = 0

def get_cut_size(cell_list, net_list):
  cut_size = 0

  first_partition = None

  for net in net_list:
    toggle = 0
    for cell in net_list[net].cell_list:
      if toggle == 0:
        first_partition = cell_list[cell].cell_partition
        toggle = 1

      else:
        if cell_list[cell].cell_partition != first_partition:
          cut_size += 1
          break

  return cut_size


def generate_cell_list(path):
  Cell_list = {}
  with open(path) as file:
    file_data = file.readlines()

    for i, cell_data in enumerate(file_data):
      cell_data = cell_data.split(" ")
      cell_id = cell_data[0]
      new_cell = Cell(cell_id)
      new_cell.cell_size = int(cell_data[1])

      Cell_list[new_cell.id] = new_cell


  return Cell_list

def generate_net_list(path):
  Net_list = {}
  net_id = 0
  i = 5
  global cell_map


  with open(path) as file:
    file_data = file.readlines()
    length = len(file_data)
    while i < length:

      net_id += 1
      new_net = Net(net_id)
      cell_id = file_data[i].split(" ")[0]
      new_net.cell_list.append(cell_id)
      cell_map[cell_id].net_list.add(net_id)
      if cell_map[cell_id].cell_partition == 0:
          new_net.Asize += 1
      else:
          new_net.Bsize += 1
      i+= 1

      while (i< length and (file_data[i].split(" ")[1] != "s")):

        cell_id = file_data[i].split(" ")[0]
        new_net.cell_list.append(cell_id)
        cell_map[cell_id].net_list.add(net_id)
        if cell_map[cell_id].cell_partition == 0:
            new_net.Asize += 1
        else:
            new_net.Bsize += 1
        i+= 1

      Net_list[net_id] = new_net


  return Net_list

def initial_partitions():
  global cell_map, sizeA, sizeB
  for node in cell_map:
    partition = random.choice([0, -1])
    if partition == 0:
        if sizeB + cell_map[node].cell_size - sizeA <= int(len(cell_map) * (0.1)):
            partition = ~partition
            sizeB += cell_map[node].cell_size
    elif partition == -1:
        if sizeA + cell_map[node].cell_size - sizeB <= int(len(cell_map) * (0.1)):
            partition = ~partition
            sizeA += cell_map[node].cell_size
    cell_map[node].cell_partition = partition
    cell_map[node].lock_status = 0


def isbalanced(target_cell):
    global sizeA, sizeB, cell_map
    temp_sizeA = sizeA
    temp_sizeB = sizeB
    if cell_map[target_cell].cell_partition == 0:
        temp_sizeA -= cell_map[target_cell].cell_size
        temp_sizeB += cell_map[target_cell].cell_size
    elif cell_map[target_cell].cell_partition == -1:
        temp_sizeA += cell_map[target_cell].cell_size
        temp_sizeB -= cell_map[target_cell].cell_size

    if abs(temp_sizeA - temp_sizeB) <= int(len(cell_map) * (0.1)):
        sizeA = temp_sizeA
        sizeB = temp_sizeB
        return True
    else:
        return False

def update_buckets(curr_cell, prev_gain):
    global gain_bucket
    new_gain = cell_map[curr_cell].cell_gain
    gain_bucket[prev_gain].remove(curr_cell)
    if not gain_bucket[prev_gain]:
        del gain_bucket[prev_gain]
    gain_bucket.setdefault(new_gain, []).append(curr_cell)

def update_gains(target_cell):
    global old_partition, gain_bucket
    for curr_net in cell_map[target_cell].net_list:
        if old_partition == 0:
            cell_map[target_cell].F = net_map[curr_net].Asize
            cell_map[target_cell].T = net_map[curr_net].Bsize
        elif old_partition == -1:
            cell_map[target_cell].F = net_map[curr_net].Bsize
            cell_map[target_cell].T = net_map[curr_net].Asize

        if cell_map[target_cell].T == 0:
            for curr_cell in net_map[curr_net].cell_list:
                if not cell_map[curr_cell].lock_status:
                    prev_gain = cell_map[curr_cell].cell_gain
                    cell_map[curr_cell].cell_gain += 1
                    update_buckets(curr_cell, prev_gain)

        elif cell_map[target_cell].T == 1:
            for curr_cell in net_map[curr_net].cell_list:
                if cell_map[target_cell].cell_partition == cell_map[curr_cell].cell_partition:
                    if not cell_map[curr_cell].lock_status:
                        prev_gain = cell_map[curr_cell].cell_gain
                        cell_map[curr_cell].cell_gain -= 1
                        update_buckets(curr_cell, prev_gain)
                        break

        cell_map[target_cell].F -= 1
        cell_map[target_cell].T += 1

        if old_partition == -1:
            net_map[curr_net].Asize = cell_map[target_cell].T
            net_map[curr_net].Bsize = cell_map[target_cell].F
        elif old_partition == 0:
            net_map[curr_net].Bsize = cell_map[target_cell].T
            net_map[curr_net].Asize = cell_map[target_cell].F

def fiducciaMathAlgo():
    global best_cut, sizeA, sizeB, cutset, minimum_cut, startcut, old_partition, locked_cells, gain_bucket, max_gain_index, sizeA_minimum_cut, sizeB_minimum_cut, locked_cell_index

    passes = 1

    cutset = get_cut_size(cell_map, net_map)
    startcut = cutset
    minimum_cut = cutset
    best_cut = float("inf")

    print("[INFO] INITIAL CUT SIZE:", cutset)
    

    for pass_num in range(1, passes + 1):
        if pass_num > 1:
            cutset = minimum_cut
            sizeA = sizeA_minimum_cut
            sizeB = sizeB_minimum_cut
            for i in range(len(locked_cells)):
                curr_cell = locked_cells[i]
                cell_map[curr_cell].lock_status = 0
                if locked_cell_index <= i:
                    cell_map[curr_cell].cell_partition = ~cell_map[curr_cell].cell_partition

        gain_bucket.clear()
        for curr_cell in cell_map:
            cell_map[curr_cell].cell_gain = 0
            for curr_net in cell_map[curr_cell].net_list:
                if cell_map[curr_cell].cell_partition == 0:
                    cell_map[curr_cell].F = net_map[curr_net].Asize
                    cell_map[curr_cell].T = net_map[curr_net].Bsize
                elif cell_map[curr_cell].cell_partition == -1:
                    cell_map[curr_cell].F = net_map[curr_net].Bsize
                    cell_map[curr_cell].T = net_map[curr_net].Asize
                if cell_map[curr_cell].F == 1:
                    cell_map[curr_cell].cell_gain += 1
                if cell_map[curr_cell].T == 0:
                    cell_map[curr_cell].cell_gain -= 1
            gain_bucket.setdefault(cell_map[curr_cell].cell_gain, []).append(curr_cell)

        locked_cells.clear()
        max_gain_index = max(gain_bucket.keys())

        while True:
            toggle = True
            if not gain_bucket:
                break
            else:
                target_cell = gain_bucket[max_gain_index][-1]
                temp_index = 0
            while toggle and not isbalanced(target_cell):
                if target_cell == gain_bucket[max_gain_index][0]:
                    if max_gain_index == min(gain_bucket):
                        break
                    max_gain_index = max_gain_index - 1
                    while max_gain_index > min(gain_bucket.keys()) and gain_bucket.get(max_gain_index) is None:
                        max_gain_index -= 1
                    temp_index = 0
                    toggle = False
                    continue
                temp_index = temp_index + 1
                target_cell = gain_bucket[max_gain_index][-1 - temp_index]

            if not toggle:
              continue

            cell_map[target_cell].lock_status = 1
            locked_cells.append(target_cell)
            new_cutset = cutset - cell_map[target_cell].cell_gain
            gain_bucket[max_gain_index].remove(target_cell)
            if not gain_bucket[max_gain_index]:
                del gain_bucket[max_gain_index]
            cutset = new_cutset
            old_partition = cell_map[target_cell].cell_partition
            cell_map[target_cell].cell_partition = ~cell_map[target_cell].cell_partition
            update_gains(target_cell)
            if 0 < new_cutset <= minimum_cut:
                minimum_cut = new_cutset
                sizeA_minimum_cut = sizeA
                sizeB_minimum_cut = sizeB
                locked_cell_index = len(locked_cells)
            if gain_bucket:
                max_gain_index = max(gain_bucket.keys())
        # current_cut = get_cut_size(cell_map, net_map)
        if minimum_cut < best_cut:
          best_cut = minimum_cut
        else:
          break
    return best_cut    

if __name__ == "__main__":

    datasets = [
        "ibm01", "ibm02", "ibm03", "ibm04", "ibm05", "ibm06",
        "ibm07", "ibm08", "ibm09", "ibm10", "ibm11", "ibm12",
        "ibm13", "ibm14", "ibm15", "ibm16", "ibm17", "ibm18",
    ]
    current_directory = os.getcwd()

    for ibm_file in datasets:
        print("[INFO] Current Benchmark:", ibm_file)

        startime = datetime.datetime.now()
        print("[INFO] STARTTIME:", startime)

        cell_map = generate_cell_list(current_directory+'/dataset/'+ibm_file+'.are')
        initial_partitions()
        net_map = generate_net_list(current_directory+'/dataset/'+ibm_file+'.net')

        best_cut = fiducciaMathAlgo()
        print("[INFO] BEST CUT:", best_cut)

        endtime = datetime.datetime.now()
        print("[INFO] ENDTIME : ", endtime)

        print("[INFO] Total Execution: ", endtime - startime)

        print("\n")
