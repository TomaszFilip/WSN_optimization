from math import sqrt
import random
import networkx as nx
import matplotlib.pyplot as plt
import sys

class Sensor:
    def __init__(self, id, e_max, p, x_coordination, y_coordination, z_coordination):
        self.id = id
        self.energy = e_max
        self.p = p
        self.Ti = 0
        self.x_coordination = x_coordination
        self.y_coordination = y_coordination
        self.z_coordination = z_coordination
        self.L = None
        self.M = None
        self.previous = None
    
    def show(self):
        print(f"==============\nSensor\nId:{self.id}\nx: {self.x_coordination}, y: {self.y_coordination}, z: {self.z_coordination}")


class Link:
    def __init__(self, first_sensor: Sensor, second_sensor : Sensor, distance):
        self.first_sensor = first_sensor
        self.second_sensor = second_sensor
        self.distance = distance
    
    def show(self):
        print(f"==============\nLink\nFirstSensorId:{self.first_sensor.id}\nSecondSensorId: {self.second_sensor.id}\nDistance: {self.distance}")


class Graph:

    def __init__(self, sensors, links):
        self.sensors = sensors
        self.links = links
        

def create_graph(sensors_number: int, e_max: float, p: float, min_x: float, 
                 max_x: float, min_y: float, max_y: float, min_z: float, max_z: float, visualize: bool):
    sensors = []
    links  = []
    for i in range(sensors_number):
        x = random.randint(min_x, max_x)
        y = random.randint(min_y, max_y)
        z = random.randint(min_z, max_z)

        sensors.append(Sensor(id= i + 1, e_max = e_max, p = p, x_coordination = x, y_coordination = y, 
                              z_coordination = z))
    
    return_links = []
    for i in range(1, sensors_number + 1):
        for j in range(i + 1, sensors_number + 1):
            distance  = sqrt(pow(sensors[i - 1].x_coordination - sensors[j - 1].x_coordination, 2) + 
                             pow(sensors[i - 1].y_coordination - sensors[j - 1].y_coordination, 2) + 
                             pow(sensors[i - 1].z_coordination - sensors[j - 1].z_coordination, 2))
            links.append(Link(first_sensor = sensors[i - 1], second_sensor = sensors[j - 1], distance = distance))
            return_links.append(Link(second_sensor = sensors[i - 1], first_sensor = sensors[j - 1], distance = distance))
    links = links + return_links
            
    #Visualize graph

    if visualize is True:

        for sensor in sensors: 
            sensor.show()

        for link in links[:int(len(links) / 2)]:
            link.show()

        nodes_v = [str(x) for x in range(1, sensors_number + 1)] 
        edges_v = [(str(x.first_sensor.id), str(x.second_sensor.id)) for x in links[:int(len(links) / 2)]]
        edges_labels_v = {(str(x.first_sensor.id), str(x.second_sensor.id)): str(round(x.distance, 2)) for x in links[:int(len(links) / 2)]}

        G = nx.Graph()
        G.add_nodes_from(nodes_v)
        G.add_edges_from(edges_v)

        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(G)
        nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=500, font_size=10, 
                font_weight='bold', edge_color='gray')

        nx.draw_networkx_edge_labels(G, pos, edge_labels=edges_labels_v, font_color='red', font_size=9)

        plt.title("Network graph")
        plt.show()
    
    return Graph(sensors, links)

def neighbour(G, W, r,e_min):
    W_neighbours = list(set(set(G.sensors) - set(W)))

    available_neighbours = []

    for link in G.links:
        if link.first_sensor in W and link.second_sensor in W_neighbours:
            if r >= link.distance and link.second_sensor.energy >= e_min:
                available_neighbours.append(link.second_sensor)

    return available_neighbours

def dijkstra(graph :Graph, vs, vq, r, e_min):
    #Initializing

    W = []  #W = Null
    graph.sensors[vq - 1].L = 0

    for x in range(len(graph.sensors)):
        graph.sensors[x].L = sys.maxsize 
        graph.sensors[x].M = 1
   
      
    W.append(graph.sensors[vq - 1])
    
    while vs not in [w.id for w in W]  and len(neighbour(graph, W, r, e_min)) > 0:
        


        pass
        




        
        
   # return

def utility_function(node: Sensor, next_node:Sensor, starting_node, end_node, beta, gamma, T, Q, h, Emax, product):
    C=calc_Cij(node, next_node, beta, gamma, starting_node, T,Q, h, Emax)
    return product-C

def get_distance(node:Sensor, nextNode:Sensor, links :set[Link]):
    for link in links:
        if (link.first_sensor.id == node.id and link.second_sensor.id == nextNode.id):
            return link.distance
    return None  


#Q - transmission payment (equal for all transmissions)
#T - number of all transmissions in the network
#h - number of nodes on path
#Emax
#M - calculated based on probabilities on the path
def calc_Cij(node:Sensor, nextNode:Sensor, beta :float, gamma:float, starting_node : int, T:int, Q:float, h:int, Emax:float):
    Dij=get_distance(node,nextNode)

    if node.id==starting_node:
        Cij=beta/(node.M-h*Q) * Dij*Dij * Emax/nextNode.energy * (1+gamma*nextNode.Ti/T)
    else:
        Cij=beta/Q * Dij*Dij * Emax/nextNode.energy * (1+gamma*nextNode.Ti/T)
    return Cij


# random_parameters
min_x = 30
max_x = 70
min_y = 30
max_y = 70
min_z = 30
max_z = 70
# visualization
visualize = False

#algorithm parameters
sensors_number = 5
e_max = 100
e_min = 30
p = 85
r = 50
starting_node = 1
end_node = 3

graph = create_graph(sensors_number, e_max, p, min_x, max_x, min_y, max_y, min_z, max_z, visualize)
dijkstra(graph, starting_node, end_node, r, e_min)