#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>
using namespace std;

const int WIDTH = 500; // Width of the canvas
const int HEIGHT = 500; // Height of the canvas
const double VIEWPORT_W = 1; // Width of the viewport
const double VIEWPORT_H = 2; // Height of the viewport
const double PROJECTION_PLANE_Z = 1; // Z-coordinate of the projection plane
const int MAX_RECURSION_DEPTH = 3; // Maximum recursion depth for reflections

// Structs for 3D vectors, colors, spheres, planes, and lights
struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(double s) const { return Vec3(x / s, y / s, z / s); }
    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 normalize() const {
        double len = sqrt(x * x + y * y + z * z);
        return (len == 0) ? *this : *this / len;
    }
};

struct Color {
    int r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(int r, int g, int b) : r(r), g(g), b(b) {}
    Color operator*(double s) const { return Color(min(255, int(r * s)), min(255, int(g * s)), min(255, int(b * s))); }
    Color operator+(const Color& c) const { return Color(min(255, r + c.r), min(255, g + c.g), min(255, b + c.b)); }
};

struct Sphere {
    Vec3 center;
    double radius;
    Color color;
    int specular;
    double reflective;
};

struct Plane {
    Vec3 point;
    Vec3 normal;
    Color color;
    int specular;
    double reflective;
};

struct Light {
    string type;
    double intensity;
    Vec3 position;
    Vec3 direction;
};

// Global variables for spheres, planes, lights, camera position, and background color
vector<Sphere> spheres;
vector<Plane> planes;
vector<Light> lights;
Vec3 camera(0, 0, -1.0);
Color BACKGROUND(0, 0, 0);

Vec3 CanvasToViewport(int x, int y) { // Convert canvas coordinates to viewport coordinates
    return Vec3(x * VIEWPORT_W / WIDTH, y * VIEWPORT_H / HEIGHT, PROJECTION_PLANE_Z);
}

// Function to check if a ray intersects with a sphere
pair<double, double> IntersectRaySphere(const Vec3& O, const Vec3& D, const Sphere& sphere) {
    Vec3 CO = O - sphere.center; // Center of the sphere
    double a = D.dot(D); // Direction of the ray
    double b = 2 * CO.dot(D); // Coefficient for the quadratic equation
    double c = CO.dot(CO) - sphere.radius * sphere.radius; // Distance from the ray to the center of the sphere
    double discriminant = b * b - 4 * a * c; // Discriminant of the quadratic equation
    if (discriminant < 0) return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()}; // No intersection
    double t1 = (-b + sqrt(discriminant)) / (2 * a); // First intersection point
    double t2 = (-b - sqrt(discriminant)) / (2 * a); // Second intersection point
    return {t1, t2}; // Return the intersection points
}

// Function to check if a ray intersects with a plane
double IntersectRayPlane(const Vec3& O, const Vec3& D, const Plane& plane) {
    double denom = plane.normal.dot(D); // Denominator for the intersection equation
    if (abs(denom) < 1e-6) return numeric_limits<double>::infinity(); // No intersection
    double t = (plane.point - O).dot(plane.normal) / denom; // Intersection point
    return t >= 0 ? t : numeric_limits<double>::infinity(); // Return the intersection point if it's in front of the ray
}

// Struct to hold hit information
struct HitInfo {
    double t; // Distance to the intersection point
    Color color;
    Vec3 normal; // Normal at the intersection point
    int specular; // Specular exponent for the material
    double reflective; // Reflective property of the material
};

// Function to find the closest intersection of a ray with all objects in the scene
bool ClosestIntersection(const Vec3& O, const Vec3& D, double t_min, double t_max, HitInfo& hit) {
    hit.t = numeric_limits<double>::infinity(); // Initialize hit distance to infinity
    bool found = false; // Flag to check if any intersection was found

    for (const Sphere& s : spheres) { // Iterate through all spheres
        // Check for intersection with the sphere
        auto [t1, t2] = IntersectRaySphere(O, D, s); // Get intersection points with the sphere
        if (t1 >= t_min && t1 <= t_max && t1 < hit.t) { // Check if the first intersection point is valid
            hit.t = t1; // Update hit distance
            Vec3 P = O + D * t1; // Calculate intersection point
            hit.normal = (P - s.center).normalize(); // Calculate normal at the intersection point
            hit.color = s.color; // Set color of the hit object
            hit.specular = s.specular; // Set specular exponent
            hit.reflective = s.reflective; // Set reflective property
            found = true; // Mark that an intersection was found
        }
        if (t2 >= t_min && t2 <= t_max && t2 < hit.t) { // Check if the second intersection point is valid
            hit.t = t2; // Update hit distance
            Vec3 P = O + D * t2; // Calculate intersection point
            hit.normal = (P - s.center).normalize(); // Calculate normal at the intersection point
            hit.color = s.color; // Set color of the hit object
            hit.specular = s.specular; // Set specular exponent
            hit.reflective = s.reflective; // Set reflective property
            found = true; // Mark that an intersection was found
        }
    }
 // Iterate through all planes
    for (const Plane& p : planes) { // Check for intersection with the plane
        double t = IntersectRayPlane(O, D, p); // Get intersection point with the plane
        if (t >= t_min && t <= t_max && t < hit.t) { // Check if the intersection point is valid
            hit.t = t; // Update hit distance
            hit.normal = p.normal; // Set normal of the plane
            hit.color = p.color; // Set color of the hit object
            hit.specular = p.specular; // Set specular exponent
            hit.reflective = p.reflective; // Set reflective property
            found = true; // Mark that an intersection was found
        }
    }
    // Return true if any intersection was found, false otherwise
    return found;
}

Vec3 ReflectRay(const Vec3& R, const Vec3& N) { // Reflect a ray R around a normal N
    // R is the incoming ray, N is the normal at the intersection point
    return N * (2 * N.dot(R)) - R; // Reflect the ray using the normal
}

// Function to compute lighting at a point P with normal N and view direction V
double ComputeLighting(const Vec3& P, const Vec3& N, const Vec3& V, int specular) { // Compute lighting at a point P
    // N is the normal at the point, V is the view direction, specular is the specular exponent
    double intensity = 0.0; // Initialize intensity to 0
    // Iterate through all lights in the scene
    for (const Light& light : lights) {
        Vec3 L;
        double t_max;

        if (light.type == "ambient") { // Ambient light
            intensity += light.intensity;
            continue;
        } else if (light.type == "point") { // Point light
            L = light.position - P;
            t_max = 1.0;
        } else {
            L = light.direction; // Directional light
            t_max = numeric_limits<double>::infinity();
        }

        // Shadow check
        HitInfo shadowHit; // Initialize shadow hit information
        // Check if the ray from the point P to the light source intersects with any object in the scene
        if (ClosestIntersection(P, L.normalize(), 0.001, t_max, shadowHit)) continue;

        // Diffuse
        double n_dot_l = N.dot(L); // Calculate the dot product of the normal and light direction
        // If the dot product is positive, the light is hitting the surface
        if (n_dot_l > 0)
            intensity += light.intensity * n_dot_l / (N.normalize().dot(N.normalize()) * L.dot(L));

        // Specular
        if (specular != -1) {
            Vec3 R = ReflectRay(L, N); // Reflect the light direction around the normal
            // Calculate the dot product of the reflected ray and the view direction
            double r_dot_v = R.dot(V);
            if (r_dot_v > 0)
                intensity += light.intensity * pow(r_dot_v / (R.dot(R) * V.dot(V)), specular); // Calculate specular intensity
        }
    }
    return intensity;
}

Color TraceRay(const Vec3& O, const Vec3& D, double t_min, double t_max, int depth) { // Trace a ray O in direction D
    // t_min is the minimum distance to check for intersections, t_max is the maximum distance
    HitInfo hit;
    if (!ClosestIntersection(O, D, t_min, t_max, hit)) return BACKGROUND; // Check for intersections with all objects in the scene
    // If no intersection is found, return the background color

    // Calculate the intersection point and view direction
    // P is the intersection point, V is the view direction (opposite of D)
    Vec3 P = O + D * hit.t; // Calculate the intersection point
    // V is the view direction (opposite of D)
    Vec3 V = D * -1; // Calculate the view direction
    // Compute lighting at the intersection point
    double lighting = ComputeLighting(P, hit.normal, V, hit.specular); // Calculate lighting at the intersection point
    // Calculate the color at the intersection point
    Color localColor = hit.color * lighting; // Calculate the local color based on the hit object's color and lighting
    // If the depth is 0 or the object is not reflective, return the local color

    if (depth <= 0 || hit.reflective <= 0) return localColor; // If the maximum recursion depth is reached or the object is not reflective, return the local color

    // Calculate the reflected ray direction and origin
    // R is the reflected ray direction, reflect_origin is the origin of the reflected ray
    Vec3 R = ReflectRay(V, hit.normal);
    Vec3 reflect_origin = P + hit.normal * 1e-4; // Offset the origin slightly to avoid self-intersection
    // Trace the reflected ray and calculate the color at the intersection point
    Color reflectedColor = TraceRay(reflect_origin, R, 0.001, numeric_limits<double>::infinity(), depth - 1); // Trace the reflected ray recursively
    // Combine the local color and reflected color based on the reflective property of the object
    return localColor * (1 - hit.reflective) + reflectedColor * hit.reflective;
}

int main() {
    const int NUM_FRAMES = 45; // Number of frames to render
    double max_distance = 3.0; // Maximum distance for the spheres to move
    double z_offset = 0.05; // Offset for the spheres in the z-direction

    // Loop through each frame and render the scene
    for (int frame = 0; frame < NUM_FRAMES; ++frame) {
        string filename = "frame" + to_string(frame) + ".ppm"; // Output filename for the frame
        ofstream out(filename); // Open the output file for writing
        out << "P3\n" << WIDTH << " " << HEIGHT << "\n255\n"; // Write the PPM header

        // Clear the scene for each frame
        spheres.clear();
        planes.clear();
        lights.clear();

        // Define the camera position and background color
        lights.push_back({"ambient", 0.5});
        lights.push_back({"point", 0.8, Vec3(2, 1, 0)});
        lights.push_back({"directional", 0.2, Vec3(), Vec3(1, 4, 4).normalize()});

        // Define the spheres and planes in the scene
        double t = frame / double(NUM_FRAMES - 1); // Normalize the frame number to [0, 1]
        double pingpong = t <= 0.5 ? t * 2 : (1.0 - t) * 2; // Ping-pong effect for the spheres
        double angle = pingpong * M_PI; // Calculate the angle for the spheres' movement
        double x_offset = max_distance * cos(angle); // Calculate the x-offset for the spheres' movement

        // Define the spheres in the scene
        spheres.push_back({Vec3(x_offset, 0, 3 - z_offset), 0.7, Color(65, 0, 85), 500, 0.3});
        spheres.push_back({Vec3(-x_offset, 0, 3 + z_offset), 0.7, Color(50, 60, 80), 10, 0.4});

        // Define the planes in the scene
        planes.push_back({Vec3(0, -2, 0), Vec3(0, 1, 0), Color(200, 200, 200), 1000, 0.5});
        planes.push_back({Vec3(0, 2, 0), Vec3(0, -1, 0), Color(200, 200, 200), 1000, 0.5});
        

        // Loop through each pixel in the canvas and trace rays to render the scene
        for (int y = -HEIGHT / 2; y < HEIGHT / 2; ++y) { // Loop through each row of pixels
            for (int x = -WIDTH / 2; x < WIDTH / 2; ++x) { // Loop through each column of pixels
                Vec3 D = CanvasToViewport(x, y); // Convert canvas coordinates to viewport coordinates
                Color color = TraceRay(camera, D, 1.0, numeric_limits<double>::infinity(), MAX_RECURSION_DEPTH); // Trace the ray and get the color
                out << color.r << " " << color.g << " " << color.b << " "; // Write the color to the output file
            }
            out << "\n"; // New line after each row of pixels
        }

        out.close(); // Close the output file
        // Print a message indicating that the frame has been rendered
        cout << "Rendered frame " << frame << " → " << filename << endl;
    }

    return 0;
}
