''' An implementation of "UNIQUENESS OF ORTHOGONAL CONNECT-THE-DOTS" by Joseph O'Rourke
'''
import numpy as np 
import random 
import matplotlib.pyplot as plt
from timeout import timeout
from settings import EPSILON
import networkx as nx
grid_size = 5
points_size = 50

def get_dot_dict_from_dots(dots, direction):
  x_index, y_index = 0, 1
  if direction == 'v':
    x_index, y_index = 1, 0
  dot_dict = dict()
  for dot in dots:
    if dot[x_index] not in dot_dict:
      dot_dict[dot[x_index]] = [dot[y_index]]
    else:
      dot_dict[dot[x_index]].append(dot[y_index])
  return dot_dict

def remove_repeat_from_sorted_list(a):
  for i in range(len(a)-1, 0, -1):
      if a[i] == a[i-1]:
        del a[i]
  return a 

def check_dots_parity(dot_dict):
  for dot_c in dot_dict:
    row_d = remove_repeat_from_sorted_list(list(sorted(dot_dict[dot_c])))
    if len(row_d)%2 == 1:
      return False, []
    dot_dict[dot_c] = row_d
  return True, dot_dict

def check_two_edge_cross(h_edge, v_edge):
  (h_x, h_y1, h_y2) = h_edge
  (v_y, v_x1, v_x2) = v_edge
  if (h_x > v_x1) and (h_x < v_x2) and (v_y > h_y1) and (v_y < h_y2):
    return True
  return False

def generate_points(grid_size, points_size):
  # ensure the number of vertices on each line is even
  x_coords = [random.randint(1, grid_size) for x in range(points_size)]
  y_coords = [random.randint(1, grid_size) for x in range(points_size)]
  dots = zip(x_coords, y_coords)
  h_dot_dict = get_dot_dict_from_dots(dots, 'h')
  is_good, h_dot_dict = check_dots_parity(h_dot_dict)
  if not is_good:
    return [], [], []
  recover_dots = []
  for dot_x in h_dot_dict:
    for dot_y in h_dot_dict[dot_x]:
      recover_dots.append([dot_x, dot_y])
  recover_dots = np.array(recover_dots)
  v_dot_dict = get_dot_dict_from_dots(recover_dots, 'v')
  is_good, v_dot_dict = check_dots_parity(v_dot_dict)
  if not is_good:
    return [], [], []
  return recover_dots, h_dot_dict, v_dot_dict

def generate_edges(h_dot_dict, v_dot_dict):
  h_edges = []
  v_edges = []
  edges = []
  for dot_x in h_dot_dict:
    for i in range(int(len(h_dot_dict[dot_x])/2)):
      edges.extend([[dot_x, h_dot_dict[dot_x][2*i]], [dot_x, h_dot_dict[dot_x][2*i+1]], [np.nan, np.nan]])
      h_edges.append([dot_x, h_dot_dict[dot_x][2*i], h_dot_dict[dot_x][2*i+1]])
  for dot_y in v_dot_dict:
    for i in range(int(len(v_dot_dict[dot_y])/2)):
      edges.extend([[v_dot_dict[dot_y][2*i], dot_y], [v_dot_dict[dot_y][2*i+1], dot_y], [np.nan, np.nan]])
      v_edges.append([dot_y, v_dot_dict[dot_y][2*i], v_dot_dict[dot_y][2*i+1]])
  for h_edge in h_edges:
    for v_edge in v_edges:
      if check_two_edge_cross(h_edge, v_edge):
        print('cross edge! Retry...')
        return []
  edges = np.array(edges)
  return edges

@timeout(5)
def generate_orthogonal_map():
  edges = []
  while len(edges) == 0:
    dots, horizontal_dots, vertical_dots = [], [], []
    while len(dots) == 0:
      dots, horizontal_dots, vertical_dots = generate_points(grid_size, points_size)
    edges = generate_edges(horizontal_dots, vertical_dots)
  return dots, edges

def get_largest_polygon(edges):
  graph_edges = [(tuple(edges[3*i]), tuple(edges[3*i+1])) for i in range(int(len(edges)/3))]
  polygon = nx.Graph()
  polygon.add_edges_from(graph_edges)
  cycles = nx.cycle_basis(polygon)
  cycle_len = [len(cycle) for cycle in cycles]
  longest_cycle = cycles[cycle_len.index(max(cycle_len))]
  largest_polygon = [[point[0], point[1]] for point in longest_cycle]
  largest_polygon.append([longest_cycle[0][0], longest_cycle[0][1]])
  largest_polygon = np.array(largest_polygon)
  return largest_polygon, longest_cycle

def get_line_number_in_file(fname):
  count = 0
  with open(fname) as f:
    for line in f:
      if line.strip():
        count += 1
  return count

if __name__ == '__main__':
  env_file_name = './orthogonal_env.py'
  start_line = get_line_number_in_file(env_file_name)
  env_count = 0
  while env_count < 10:
    dots = []
    edges = []
    count = 0
    while(len(dots) == 0 and count < 10):
      count += 1
      dots, edges = generate_orthogonal_map()
    largest_polygon, longest_cycle = get_largest_polygon(edges)
    if len(longest_cycle) > 10:
      with open(env_file_name, 'a') as f:
        f.write('\northogonal_env_{} = np.array({})\n'.format(start_line+env_count, longest_cycle))
        env_count += 1

      plt.figure()
      ax = plt.gca()
      ax.axis('equal')
      ax.axis('off')
      plt.plot(edges[:, 0], edges[:, 1], 'b-')
      plt.plot(dots[:, 0], dots[:, 1], 'ko')
      plt.show()

      plt.figure()
      ax = plt.gca()
      ax.axis('equal')
      ax.axis('off')
      plt.plot(largest_polygon[:, 0], largest_polygon[:, 1], 'b-')
      plt.show()
