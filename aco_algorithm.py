import math
import random
import matplotlib.pyplot as plt
# specify which coordinates you want to use by providing the name of the file and variable
from coordinates_cluster_Γ import coordinates

def distance(coord1, coord2):
    lat1, lon1 = coord1[1], coord1[2]
    lat2, lon2 = coord2[1], coord2[2]
    radius = 6371  # Earth's radius in kilometers
    if coord1 == coord2:
        return 0  # Return 0 distance for the same coordinates
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) * math.sin(dlat / 2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon / 2) * math.sin(dlon / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = radius * c
    return distance

# ACO algorithm implementation
class ACO:
    def __init__(self, coordinates, alpha=1.0, beta=5.0, evaporation_rate=0.5, num_ants=100, num_iterations=300):
        self.coordinates = coordinates
        self.num_cities = len(coordinates)
        self.distances = [[distance(coord1, coord2) + 1e-10 for coord2 in coordinates] for coord1 in coordinates]
        self.pheromones = [[1 / distance for distance in distances] for distances in self.distances]
        self.alpha = alpha
        self.beta = beta
        self.evaporation_rate = evaporation_rate
        self.num_ants = num_ants
        self.num_iterations = num_iterations

    def find_optimal_path(self, start_city, end_city):
        optimal_path = []
        optimal_distance = float('inf')
        for _ in range(self.num_iterations):
            paths = self.construct_paths(start_city, end_city)
            path_distances = [self.calculate_path_distance(path) for path in paths]
            min_distance = min(path_distances)
            if min_distance < optimal_distance:
                optimal_distance = min_distance
                optimal_path = paths[path_distances.index(min_distance)]
            self.update_pheromones(paths, path_distances)
        return optimal_path, optimal_distance

    def construct_paths(self, start_city, end_city):
        paths = []
        for _ in range(self.num_ants):
            unvisited = set(range(self.num_cities))
            unvisited.remove(start_city)
            unvisited.remove(end_city)
            current_city = start_city
            path = [current_city]
            while unvisited:
                next_city = self.construct_path(current_city, unvisited)
                path.append(next_city)
                unvisited.remove(next_city)
                current_city = next_city
            path.append(end_city)
            paths.append(path)
        return paths

    def construct_path(self, current_city, unvisited):
        total_probabilities = 0
        probabilities = []
        for city in unvisited:
            pheromone = self.pheromones[current_city][city]
            distance = self.distances[current_city][city]
            probability = math.pow(pheromone, self.alpha) * math.pow(1 / distance, self.beta)
            probabilities.append(probability)
            total_probabilities += probability
        if total_probabilities == 0:
            return random.choice(list(unvisited))
        probabilities = [probability / total_probabilities for probability in probabilities]
        next_city = random.choices(list(unvisited), probabilities)[0]
        return next_city

    def calculate_path_distance(self, path):
        total_distance = 0
        for i in range(len(path) - 1):
            total_distance += self.distances[path[i]][path[i + 1]]
        return total_distance

    def update_pheromones(self, paths, path_distances):
        for i in range(self.num_cities):
            for j in range(self.num_cities):
                self.pheromones[i][j] *= self.evaporation_rate
        for path, distance in zip(paths, path_distances):
            for i in range(len(path) - 1):
                city1 = path[i]
                city2 = path[i + 1]
                self.pheromones[city1][city2] += 1.0 / distance

    def plot_optimal_path(self, path):
        fig, ax = plt.subplots()
        for i in range(len(path) - 1):
            ax.plot([self.coordinates[path[i]][2], self.coordinates[path[i + 1]][2]],
                    [self.coordinates[path[i]][1], self.coordinates[path[i + 1]][1]], 'b-')
        for coord in self.coordinates:
            ax.plot(coord[2], coord[1], 'ro', markersize=8)
            ax.annotate(coord[0], (coord[2], coord[1]), fontsize=8, ha='center', va='bottom')
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Optimal Path for Γ and best entry point to Δ")
        plt.show()

# Usage example
aco = ACO(coordinates)
start_city = 12  # Index of the start city in the coordinates list
end_city = 13  # Index of the end city in the coordinates list
optimal_path, optimal_distance = aco.find_optimal_path(start_city, end_city)
print("Optimal Path:")
for city in optimal_path:
    print(coordinates[city][0])
print("Optimal Distance:", optimal_distance, "km")
aco.plot_optimal_path(optimal_path)




