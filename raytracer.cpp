#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <thread>


using namespace parser;


void print(Vec3f v3) {
    std::cout << v3.x << " " << v3.y << " " << v3.z << std::endl;
}

void print(Vec3i v3) {
    std::cout << v3.x << " " << v3.y << " " << v3.z << std::endl;
}

float det(float a, float b, float c, float d, float e, float f, float g, float h, float i) {
    return a*(e*i-h*f) + b*(g*f-d*i) + c*(d*h-e*g);
}

float dist_btw_points(Vec3f p1, Vec3f p2) {
    return sqrt( pow((p1.x-p2.x), 2) + pow((p1.y-p2.y), 2) + pow((p1.z-p2.z), 2) );
}

Ray generateRay(Camera camera, int i, int j) {
	Ray result;
    float left=camera.near_plane.x, right=camera.near_plane.y;
    float bottom=camera.near_plane.z, top=camera.near_plane.w;
    int nx = camera.image_width;
    int ny = camera.image_height;
    Vec3f gaze = camera.gaze;
    Vec3f e = camera.position;
    float dist = camera.near_distance;
    Vec3f v = camera.up, u = Cross(gaze, v), w = Cross(u, v);

	float su = (i+0.5)*(right-left)/nx;
	float sv = (j+0.5)*(top-bottom)/ny;
	
	Vec3f m = e + gaze*dist;
    Vec3f q = m + u*left + v*top;
	Vec3f s = q + u*su + v*(-sv);

	result.o = e;
	result.d = s + e*(-1);
	
	return result;
}

float intersectSphere(Ray r, Sphere &s, std::vector<Vec3f> &vertex_data) {	
	int vertex_id = s.center_vertex_id;
    Vec3f c = vertex_data[vertex_id-1];
	
	float t,t1,t2;
	
	float C = (r.o.x-c.x)*(r.o.x-c.x)+(r.o.y-c.y)*(r.o.y-c.y)+(r.o.z-c.z)*(r.o.z-c.z)-s.radius*s.radius;

	float B = 2*r.d.x*(r.o.x-c.x)+2*r.d.y*(r.o.y-c.y)+2*r.d.z*(r.o.z-c.z);
	
	float A = r.d.x*r.d.x+r.d.y*r.d.y+r.d.z*r.d.z;
	
	float delta = B*B-4*A*C;
	
	if (delta<0) return -1;
	else if (delta==0)
	{
		t = -B / (2*A);
	}
	else
	{
		delta = sqrt(delta);
		A = 2*A;
		t1 = (-B + delta) / A;
		t2 = (-B - delta) / A;
				
		if (t1<t2) t=t1; else t=t2;
	}
	
	return t;
}

float intersectTriangle(Ray r, Face &face, std::vector<Vec3f> &vertex_data) {
    int vertex_id0 = face.v0_id;
    int vertex_id1 = face.v1_id;
    int vertex_id2 = face.v2_id;
    Vec3f a = vertex_data[vertex_id0-1];
    Vec3f b = vertex_data[vertex_id1-1];
    Vec3f c = vertex_data[vertex_id2-1];

    float detA = det(a.x-b.x, a.y-b.y, a.z-b.z, a.x-c.x, a.y-c.y, a.z-c.z, r.d.x, r.d.y, r.d.z);

    float beta = det(a.x-r.o.x, a.y-r.o.y, a.z-r.o.z, a.x-c.x, a.y-c.y, a.z-c.z, r.d.x, r.d.y, r.d.z);
    beta = beta / detA;

    float gama = det(a.x-b.x, a.y-b.y, a.z-b.z, a.x-r.o.x, a.y-r.o.y, a.z-r.o.z, r.d.x, r.d.y, r.d.z);
    gama = gama / detA;

    float t = det(a.x-b.x, a.y-b.y, a.z-b.z, a.x-c.x, a.y-c.y, a.z-c.z, a.x-r.o.x, a.y-r.o.y, a.z-r.o.z);
    t = t / detA;

    if (0 <= beta && 0 <= gama && beta+gama <= 1) {
        return t;
    }
    else {
        return -1;
    }
}

Vec3i computeColor(Ray r, Scene &scene, int rec) {
	Vec3i c;
	float minT = 900000; // some large number
	float t;
    Vec3f dif_ref, amb_ref, spec_ref, mir_ref;
    int mesh_triangle_sphere = -1; // 0-1-2
    int mesh_face_idx = -1;
	
	c.x = rec==0 ? scene.background_color.x : 0;
    c.y = rec==0 ? scene.background_color.y : 0;
    c.z = rec==0 ? scene.background_color.z : 0;
	int minI = -1;

    bool is_mirror;
    
	for (int i=0; i<scene.spheres.size(); i++) { // intersect spheres
		t = intersectSphere(r, scene.spheres[i], scene.vertex_data);
		if (t<minT && t>=0)
		{
            int material_id = scene.spheres[i].material_id;
            dif_ref = scene.materials[material_id-1].diffuse;
            amb_ref = scene.materials[material_id-1].ambient;
            spec_ref = scene.materials[material_id-1].specular;
            mir_ref = scene.materials[material_id-1].mirror;
			minI = i;
			minT = t;
            mesh_triangle_sphere = 2;
            is_mirror = scene.materials[material_id-1].is_mirror;
		}
	}

    for (int i=0; i<scene.triangles.size(); i++) { // intersect triangles
        t = intersectTriangle(r, scene.triangles[i].indices, scene.vertex_data);
        if (t<minT && t>=0) {
            int material_id = scene.triangles[i].material_id;
            dif_ref = scene.materials[material_id-1].diffuse;
            amb_ref = scene.materials[material_id-1].ambient;
            spec_ref = scene.materials[material_id-1].specular;
            mir_ref = scene.materials[material_id-1].mirror;
			minI = i;
			minT = t;
            mesh_triangle_sphere = 1;
            is_mirror = scene.materials[material_id-1].is_mirror;
        }
    }

    for (int i=0; i<scene.meshes.size(); i++) { // intersect meshes
        for (int idx=0; idx<scene.meshes[i].faces.size(); idx++) {
            t = intersectTriangle(r, scene.meshes[i].faces[idx], scene.vertex_data);
            if (t<minT && t>=0) {
                int material_id = scene.meshes[i].material_id;
                dif_ref = scene.materials[material_id-1].diffuse;
                amb_ref = scene.materials[material_id-1].ambient;
                spec_ref = scene.materials[material_id-1].specular;
                mir_ref = scene.materials[material_id-1].mirror;
                minI = i;
                minT = t;
                mesh_face_idx = idx;
                mesh_triangle_sphere = 0;
                is_mirror = scene.materials[material_id-1].is_mirror;
            }
        }
    }

	if (minI != -1) // object found
	{
        Vec3f ambient_radiance = scene.ambient_light;
        
        Vec3f P = r.o + (r.d)*minT;
        
        
        Vec3f N;
        int material_index = 0;


        // normal calculations
        if (mesh_triangle_sphere == 2) { // sphere
            int vertex_id = scene.spheres[minI].center_vertex_id;
            Vec3f center = scene.vertex_data[vertex_id-1];
            N = P + (center)*(-1);
            material_index = scene.spheres[minI].material_id;

        }

        else if (mesh_triangle_sphere == 1) { // triangle
            int vertex_id0 = scene.triangles[minI].indices.v0_id;
            int vertex_id1 = scene.triangles[minI].indices.v1_id;
            int vertex_id2 = scene.triangles[minI].indices.v2_id;
            Vec3f a = scene.vertex_data[vertex_id0-1];
            Vec3f b = scene.vertex_data[vertex_id1-1];
            Vec3f c = scene.vertex_data[vertex_id2-1];

            N = Cross(c-b, a-b);

            material_index = scene.triangles[minI].material_id;

        }

        else if (mesh_triangle_sphere == 0) { // mesh
            int vertex_id0 = scene.meshes[minI].faces[mesh_face_idx].v0_id;
            int vertex_id1 = scene.meshes[minI].faces[mesh_face_idx].v1_id;
            int vertex_id2 = scene.meshes[minI].faces[mesh_face_idx].v2_id;
            Vec3f a = scene.vertex_data[vertex_id0-1];
            Vec3f b = scene.vertex_data[vertex_id1-1];
            Vec3f c = scene.vertex_data[vertex_id2-1];

            N = Cross(c-b, a-b);

            material_index = scene.meshes[minI].material_id;

            
        }

        N = Normalized(N);

        c.x = (int) fmin((amb_ref.x * ambient_radiance.x), 255); // ambient
        c.y = (int) fmin((amb_ref.y * ambient_radiance.y), 255);
        c.z = (int) fmin((amb_ref.z * ambient_radiance.z), 255);

        for (int light_ind=0; light_ind<scene.point_lights.size(); light_ind++) {
            Vec3f light_pos = scene.point_lights[light_ind].position;
            Vec3f light_intensity = scene.point_lights[light_ind].intensity;
            float light_dist = dist_btw_points(light_pos, P);
            Vec3f L = light_pos - P;
            L = Normalized(L);

            Ray lightray;
            lightray.o = P + N*scene.shadow_ray_epsilon;
            lightray.d = light_pos - P;
            float t_light = 1;


            bool light_touches = true;

            for (int i=0; i<scene.spheres.size(); i++) { // intersect light ray with spheres
                float t_temp = intersectSphere(lightray, scene.spheres[i], scene.vertex_data);
                if (t_temp < t_light && t_temp > 0)
                {
                    light_touches = false;
                    break;
                }
            }
            
            if (light_touches == true) {
                for (int i=0; i<scene.triangles.size(); i++) { // intersect light ray with triangles
                    float t_temp = intersectTriangle(lightray, scene.triangles[i].indices, scene.vertex_data);
                    if (t_temp < t_light && t_temp > 0) {
                        light_touches = false;
                        break;
                    }
                }
            }

            if (light_touches == true) {
                for (int i=0; i<scene.meshes.size(); i++) { // intersect light ray with meshes
                    for (int idx=0; idx<scene.meshes[i].faces.size(); idx++) {
                        float t_temp = intersectTriangle(lightray, scene.meshes[i].faces[idx], scene.vertex_data);
                        if (t_temp < t_light && t_temp > 0) {
                            light_touches = false;
                            break;
                        }
                    }
                }
            }

            if (light_touches) { // if the light reaches to object, apply diffuse and specular
                Vec3f temp = dif_ref * (float) fmax(Dot(L,N), 0); // diffuse
                Vec3f irradiance = light_intensity / (float) pow(light_dist, 2);
            
                temp.x = temp.x * irradiance.x;
                temp.y = temp.y * irradiance.y;
                temp.z = temp.z * irradiance.z;
        
                c.x = (int) fmin((c.x + temp.x), 255);
                c.y = (int) fmin((c.y + temp.y), 255);
                c.z = (int) fmin((c.z + temp.z), 255);

                Vec3f p_to_camera = r.o - P; //specular

                p_to_camera = Normalized(p_to_camera);
                float phong_exp = scene.materials[material_index-1].phong_exponent;
                Vec3f half = Normalized(L + p_to_camera);

                Vec3f L_spec;
                L_spec.x = 0; L_spec.y = 0; L_spec.z = 0;
                float angle_btw_L_N = acos(Dot(L, N) / ( Magnitude(L) * Magnitude(N) ));
                if (angle_btw_L_N < (M_PI/2)) {
                    L_spec = spec_ref * ( (float) pow(fmax(Dot(N, half), 0), phong_exp) );
                    L_spec.x = L_spec.x * irradiance.x;
                    L_spec.y = L_spec.y * irradiance.y;
                    L_spec.z = L_spec.z * irradiance.z;
                }

                c.x = (int) fmin(c.x + L_spec.x, 255);
                c.y = (int) fmin(c.y + L_spec.y, 255);
                c.z = (int) fmin(c.z + L_spec.z, 255);
            }
        }

        if (is_mirror && rec < scene.max_recursion_depth) { // reflectance
            Vec3f p_to_cam = r.o - P;
            Vec3f reflected = (N * (2*(Dot(N, p_to_cam)))) - p_to_cam;
            Ray ref_ray;
            ref_ray.o = P + (N*scene.shadow_ray_epsilon);
            ref_ray.d = reflected;

            Vec3i c_ref = computeColor(ref_ray, scene, rec+1);
            c.x = (int) fmin(c.x + (c_ref.x * mir_ref.x), 255);
            c.y = (int) fmin(c.y + (c_ref.y * mir_ref.y), 255);
            c.z = (int) fmin(c.z + (c_ref.z * mir_ref.z), 255);
        }

    }
	return c;
}

void render(int width_s, int height_s, int width_e, int height_e, int i, unsigned char *image, Scene scene, int cam_idx) {
    for (int y = height_s; y < height_e; ++y)
    {
        for (int x = width_s; x < width_e; ++x)
        {
            Ray r = generateRay(scene.cameras[cam_idx], x, y);
            Vec3i color = computeColor(r, scene, 0);
            
            image[i++] = color.x;
            image[i++] = color.y;
            image[i++] = color.z;
        }
    }
}

int main(int argc, char* argv[])
{   
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    for (int cam_idx=0; cam_idx<scene.cameras.size(); cam_idx++) {
        int width = scene.cameras[cam_idx].image_width;
        int height = scene.cameras[cam_idx].image_height;

        unsigned char* image = new unsigned char [width * height * 3];

        std::thread t1(render, 0, 0, width, height/4, 0, image, scene, cam_idx);
        std::thread t2(render, 0, height/4, width, height/2, width*height*3/4, image, scene, cam_idx);
        std::thread t3(render, 0, height/2, width, (3*height)/4, width*height*3/2, image, scene, cam_idx);
        std::thread t4(render, 0, (3*height)/4, width, height, 3*width*height*3/4, image, scene, cam_idx);

        t1.join();
        t2.join();
        t3.join();
        t4.join();

        std::string filename = scene.cameras[cam_idx].image_name;
        write_ppm(filename.c_str(), image, width, height);
    }
}

