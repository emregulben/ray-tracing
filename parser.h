#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f
    {
        float x, y, z;

    };

    struct Vec3i
    {
        int x, y, z;
    };

    struct Vec4f
    {
        float x, y, z, w;
    };

    struct Ray
    {
        Vec3f o, d;
    };



    Vec3f operator/(const Vec3f& v, const float& f);
    Vec4f operator/(const Vec4f& v, const float& f);
    Vec3f operator/(const Vec3i& v, const float& f);

    Vec3f operator-(const Vec3f& v, const float& f);
    Vec4f operator-(const Vec4f& v, const float& f);
    Vec3f operator-(const Vec3i& v, const float& f);

    Vec3f operator+(const Vec3f& v, const float& f);
    Vec4f operator+(const Vec4f& v, const float& f);
    Vec3f operator+(const Vec3i& v, const float& f);

    Vec3f operator*(const Vec3f& v, const float& f);
    Vec4f operator*(const Vec4f& v, const float& f);
    Vec3f operator*(const Vec3i& v, const float& f);

    Vec3f operator+(const Vec3f& lhs, const Vec3f& rhs);
    Vec4f operator+(const Vec4f& lhs, const Vec4f& rhs);
    Vec3i operator+(const Vec3i& lhs, const Vec3i& rhs);
    Vec3f operator-(const Vec3f& lhs, const Vec3f& rhs);
    Vec4f operator-(const Vec4f& lhs, const Vec4f& rhs);
    Vec3i operator-(const Vec3i& lhs, const Vec3i& rhs);

    float Dot(const Vec3f& lhs, const Vec3f& rhs);
    float Dot(const Vec4f& lhs, const Vec4f& rhs);
    float Dot(const Vec3i& lhs, const Vec3i& rhs);
    Vec3f Cross(const Vec3f& lhs, const Vec3f& rhs);
    Vec3f Cross(const Vec3i& lhs, const Vec3i& rhs);
    float Magnitude(const Vec3f& v);
    float Magnitude(const Vec4f& v);
    float Magnitude(const Vec3i& v);
    Vec3f Normalized(const Vec3f& v);
    Vec4f Normalized(const Vec4f& v);
    Vec3f Normalized(const Vec3i& v);


    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material
    {
        bool is_mirror;
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
    };

    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;
    };

    struct Triangle
    {
        int material_id;
        Face indices;
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct Scene
    {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;

        //Functions
        void loadFromXml(const std::string &filepath);
    };
}

#endif
