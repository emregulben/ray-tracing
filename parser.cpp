#include "parser.h"
#include "tinyxml2.h"
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <iostream>

struct Vec3f {
  float x, y, z;

};
struct Vec3i {
  float x, y, z;

};
struct Vec4f {
  float x, y, z, w;
};
struct Ray {
  Vec3f o, d;
};


parser::Vec3f parser::operator/(const parser::Vec3f& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x/f;
  tmp.y = v.y/f;
  tmp.z = v.z/f;
  return tmp;
}
parser::Vec4f parser::operator/(const parser::Vec4f& v, const float& f){
  parser::Vec4f tmp;
  tmp.x = v.x/f;
  tmp.y = v.y/f;
  tmp.z = v.z/f;
  tmp.w = v.w/f;
  return tmp;
}
parser::Vec3f parser::operator/(const parser::Vec3i& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x/f;
  tmp.y = v.y/f;
  tmp.z = v.z/f;
  return tmp;
}

parser::Vec3f parser::operator-(const parser::Vec3f& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x-f;
  tmp.y = v.y-f;
  tmp.z = v.z-f;
  return tmp;
}
parser::Vec4f parser::operator-(const parser::Vec4f& v, const float& f){
  parser::Vec4f tmp;
  tmp.x = v.x-f;
  tmp.y = v.y-f;
  tmp.z = v.z-f;
  tmp.w = v.w-f;
  return tmp;
}
parser::Vec3f parser::operator-(const parser::Vec3i& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x-f;
  tmp.y = v.y-f;
  tmp.z = v.z-f;
  return tmp;
}

parser::Vec3f parser::operator+(const parser::Vec3f& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x+f;
  tmp.y = v.y+f;
  tmp.z = v.z+f;
  return tmp;
}
parser::Vec4f parser::operator+(const parser::Vec4f& v, const float& f){
  parser::Vec4f tmp;
  tmp.x = v.x+f;
  tmp.y = v.y+f;
  tmp.z = v.z+f;
  tmp.w = v.w+f;
  return tmp;
}
parser::Vec3f parser::operator+(const parser::Vec3i& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x+f;
  tmp.y = v.y+f;
  tmp.z = v.z+f;
  return tmp;
}

parser::Vec3f parser::operator*(const parser::Vec3f& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x*f;
  tmp.y = v.y*f;
  tmp.z = v.z*f;
  return tmp;
}
parser::Vec4f parser::operator*(const parser::Vec4f& v, const float& f){
  parser::Vec4f tmp;
  tmp.x = v.x*f;
  tmp.y = v.y*f;
  tmp.z = v.z*f;
  tmp.w = v.w*f;
  return tmp;
}
parser::Vec3f parser::operator*(const parser::Vec3i& v, const float& f){
  parser::Vec3f tmp;
  tmp.x = v.x*f;
  tmp.y = v.y*f;
  tmp.z = v.z*f;
  return tmp;
}

parser::Vec3f parser::operator+(const parser::Vec3f& lhs, const parser::Vec3f& rhs){
  parser::Vec3f tmp;
  tmp.x = lhs.x + rhs.x;
  tmp.y = lhs.y + rhs.y;
  tmp.z = lhs.z + rhs.z;
  return tmp;
}
parser::Vec4f parser::operator+(const parser::Vec4f& lhs, const parser::Vec4f& rhs){
  parser::Vec4f tmp;
  tmp.x = lhs.x + rhs.x;
  tmp.y = lhs.y + rhs.y;
  tmp.z = lhs.z + rhs.z;
  tmp.w = lhs.w + rhs.w;
  return tmp;
}
parser::Vec3i parser::operator+(const parser::Vec3i& lhs, const parser::Vec3i& rhs){
  parser::Vec3i tmp;
  tmp.x = lhs.x + rhs.x;
  tmp.y = lhs.y + rhs.y;
  tmp.z = lhs.z + rhs.z;
  return tmp;
}

parser::Vec3f parser::operator-(const parser::Vec3f& lhs, const parser::Vec3f& rhs){
  parser::Vec3f tmp;
  tmp.x = lhs.x - rhs.x;
  tmp.y = lhs.y - rhs.y;
  tmp.z = lhs.z - rhs.z;
  return tmp;
}
parser::Vec4f parser::operator-(const parser::Vec4f& lhs, const parser::Vec4f& rhs){
  parser::Vec4f tmp;
  tmp.x = lhs.x - rhs.x;
  tmp.y = lhs.y - rhs.y;
  tmp.z = lhs.z - rhs.z;
  tmp.w = lhs.w - rhs.w;
  return tmp;
}
parser::Vec3i parser::operator-(const parser::Vec3i& lhs, const parser::Vec3i& rhs){
  parser::Vec3i tmp;
  tmp.x = lhs.x - rhs.x;
  tmp.y = lhs.y - rhs.y;
  tmp.z = lhs.z - rhs.z;
  return tmp;
}

float parser::Dot(const parser::Vec3f& lhs, const parser::Vec3f& rhs){
  return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
}
float parser::Dot(const parser::Vec4f& lhs, const parser::Vec4f& rhs){
  return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z + lhs.w * rhs.w;
}
float parser::Dot(const parser::Vec3i& lhs, const parser::Vec3i& rhs){
  return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
}
parser::Vec3f parser::Cross(const parser::Vec3f& lhs, const parser::Vec3f& rhs){
  parser::Vec3f tmp;
  tmp.x = lhs.y*rhs.z - lhs.z*rhs.y;
  tmp.y = lhs.z*rhs.x - lhs.x*rhs.z;
  tmp.z = lhs.x*rhs.y - lhs.y*rhs.x;
  return tmp;
}
parser::Vec3f parser::Cross(const parser::Vec3i& lhs, const parser::Vec3i& rhs){
  parser::Vec3f tmp;
  tmp.x = lhs.y*rhs.z - lhs.z*rhs.y;
  tmp.y = lhs.z*rhs.x - lhs.x*rhs.z;
  tmp.z = lhs.x*rhs.y - lhs.y*rhs.x;
  return tmp;
}
    
float parser::Magnitude(const parser::Vec3f& v){
  return sqrt(pow(v.x,2) + pow(v.y, 2) + pow(v.z,2));
}

float parser::Magnitude(const parser::Vec4f& v){
  return sqrt(pow(v.x,2) + pow(v.y, 2) + pow(v.z,2) + pow(v.w, 2));
}

float parser::Magnitude(const parser::Vec3i& v){
  return sqrt(pow(v.x,2) + pow(v.y, 2) + pow(v.z,2));
}

parser::Vec3f parser::Normalized(const parser::Vec3f& v){
  return v / Magnitude(v);
}
parser::Vec4f parser::Normalized(const parser::Vec4f& v){
  return v / Magnitude(v);
}
parser::Vec3f parser::Normalized(const parser::Vec3i& v){
  return v / Magnitude(v);
}


void parser::Scene::loadFromXml(const std::string &filepath)
{
    tinyxml2::XMLDocument file;
    std::stringstream stream;

    auto res = file.LoadFile(filepath.c_str());
    if (res)
    {
        throw std::runtime_error("Error: The xml file cannot be loaded.");
    }

    auto root = file.FirstChild();
    if (!root)
    {
        throw std::runtime_error("Error: Root is not found.");
    }

    //Get BackgroundColor
    auto element = root->FirstChildElement("BackgroundColor");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0 0 0" << std::endl;
    }
    stream >> background_color.x >> background_color.y >> background_color.z;

    //Get ShadowRayEpsilon
    element = root->FirstChildElement("ShadowRayEpsilon");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0.001" << std::endl;
    }
    stream >> shadow_ray_epsilon;

    //Get MaxRecursionDepth
    element = root->FirstChildElement("MaxRecursionDepth");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0" << std::endl;
    }
    stream >> max_recursion_depth;

    //Get Cameras
    element = root->FirstChildElement("Cameras");
    element = element->FirstChildElement("Camera");
    Camera camera;
    while (element)
    {
        auto child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Gaze");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Up");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearPlane");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearDistance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageResolution");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageName");
        stream << child->GetText() << std::endl;

        stream >> camera.position.x >> camera.position.y >> camera.position.z;
        stream >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;
        stream >> camera.up.x >> camera.up.y >> camera.up.z;
        stream >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >> camera.near_plane.w;
        stream >> camera.near_distance;
        stream >> camera.image_width >> camera.image_height;
        stream >> camera.image_name;

        cameras.push_back(camera);
        element = element->NextSiblingElement("Camera");
    }

    //Get Lights
    element = root->FirstChildElement("Lights");
    auto child = element->FirstChildElement("AmbientLight");
    stream << child->GetText() << std::endl;
    stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
    element = element->FirstChildElement("PointLight");
    PointLight point_light;
    while (element)
    {
        child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;

        stream >> point_light.position.x >> point_light.position.y >> point_light.position.z;
        stream >> point_light.intensity.x >> point_light.intensity.y >> point_light.intensity.z;

        point_lights.push_back(point_light);
        element = element->NextSiblingElement("PointLight");
    }

    //Get Materials
    element = root->FirstChildElement("Materials");
    element = element->FirstChildElement("Material");
    Material material;
    while (element)
    {
        material.is_mirror = (element->Attribute("type", "mirror") != NULL);

        child = element->FirstChildElement("AmbientReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("DiffuseReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("SpecularReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("MirrorReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("PhongExponent");
        stream << child->GetText() << std::endl;

        stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
        stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
        stream >> material.specular.x >> material.specular.y >> material.specular.z;
        stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
        stream >> material.phong_exponent;

        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }

    //Get VertexData
    element = root->FirstChildElement("VertexData");
    stream << element->GetText() << std::endl;
    parser::Vec3f vertex;
    while (!(stream >> vertex.x).eof())
    {
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(vertex);
    }
    stream.clear();

    //Get Meshes
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    Mesh mesh;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> mesh.material_id;

        child = element->FirstChildElement("Faces");
        stream << child->GetText() << std::endl;
        Face face;
        while (!(stream >> face.v0_id).eof())
        {
            stream >> face.v1_id >> face.v2_id;
            mesh.faces.push_back(face);
        }
        stream.clear();

        meshes.push_back(mesh);
        mesh.faces.clear();
        element = element->NextSiblingElement("Mesh");
    }
    stream.clear();

    //Get Triangles
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    Triangle triangle;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> triangle.material_id;

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;

        triangles.push_back(triangle);
        element = element->NextSiblingElement("Triangle");
    }

    //Get Spheres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    Sphere sphere;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphere.material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphere.center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere.radius;

        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }
}


