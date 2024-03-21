//
// Changes to be made : 
// - change pow(,2) to multiply with itself
// - use all RGB channels instead of only R
// 
// 
//







/*

#include <Magick++.h> 
#include <iostream> 
#include <fstream>
#include <ctime>

using namespace std;
using namespace Magick;


void append(string txtfile, string text) {

    // reading the exising contents in the file
    fstream read(txtfile);

    // store the file contents in a variable
    string storage, buffer;
    while (getline(read, buffer)) {
        storage += buffer;
        storage += '\n';
    }

    // erase the last '\n'
    storage.erase(storage.length() - 1, 1);

    // add the new text at the end of the file
    storage += text;
    read.close();

    fstream write(txtfile);

    // write the existing contents + the added text to the file
    write << storage;
    write.close();
};


//void getPixelRGB(Image img, int x, int y, int* R, int* G, int* B) {
void getPixelRGB(Image img, int x, int y, Quantum* R, Quantum* G, Quantum* B) {
    Pixels pixel(img);

    Quantum* pix = pixel.get(x, y, 1, 1);

    *R = *pix;
    *G = *(pix + 1);
    *B = *(pix + 2);
};










int main(int argc, char** argv)
{
    InitializeMagick(*argv);

    /*
    // create a 16x16 white image with a red dot
    Image image(Geometry(16,16), Color("#007d7d"));         // #7df9ff
    Pixels p(image);
    Quantum* q = p.get(3,3,1,1);

    *q = 65535;
    *(q + 1) = 0;
    *(q + 2) = 0;

    p.sync();
    image.write("C:/Users/ACER/Desktop/asd.png");
    *






    Image image;
    image.read("C:/Users/ACER/Downloads/Picture/Batu Ferringhi.jpg");

    image.write("C:/Users/ACER/Desktop/Original.png");


    size_t width = image.size().width();
    size_t height = image.size().height();

    cout << "width   " << width << '\t' << "height   " << height << endl;

    Image resized = image;
    resized.resize(Geometry(width / 4, height / 4));
    resized.write("C:/Users/ACER/Desktop/resized.png");

    Image sampled = image;
    sampled.sample(Geometry(width / 10, height / 10));
    sampled.write("C:/Users/ACER/Desktop/sampled.png");

    Image scaled = image;
    scaled.scale(Geometry(width / 4, height / 4));
    scaled.write("C:/Users/ACER/Desktop/scaled.png");


    Image lowres;
    lowres.read("C:/Users/ACER/Desktop/sampled.png");
    lowres.sample(Geometry(width, height));
    lowres.write("C:/Users/ACER/Desktop/lowres.png");


    //cout << "After reading and writing\n";

    Image re;
    re.read("C:/Users/ACER/Desktop/lowres.png");
    size_t width_r = re.size().width();
    size_t height_r = re.size().height();



    Quantum r1, r2, g1, g2, b1, b2;
    Quantum* R1 = &r1;
    Quantum* R2 = &r2;
    Quantum* G1 = &g1;
    Quantum* G2 = &g2;
    Quantum* B1 = &b1;
    Quantum* B2 = &b2;


    if (height == height_r && width == width_r) {

        long long int cost_R = 0;

        int I = -3;
        Pixels pixel(image);
        Pixels pixel2(re);

        Quantum* pix = pixel.get(0, 0, width, height);
        Quantum* pix2 = pixel2.get(0, 0, width, height);

        for (int i = 0; i <= height - 1; i++) {
            for (int j = 0; j <= width - 1; j++) {

                I += 3;
                cost_R +=   pow((*(pix + I) - *(pix2 + I)), 2);           // multiple y it self
                //cost_G += pow((*(pix + I + 1) - *(pix2 + I + 1)), 2);
            }
            cout << i << '\t' << cost_R << endl;
        }
        cout << cost_R << endl;
    }
    else {
        cout << "Image size not compatible, not able to calculate cost function" << endl;
    }




    // lowres(/4) w original : 1434419671   (overload)
    //                         1434419671
    //                         9223372036854775807
    //                                      147713366203717
    // lowres(/10) w original long long int:147713366203717
    // lowres(/4) w original long long int: 73007418025640
    // lowres(/3) w original long long int: 44272891723260
    // lowres(/2) w original long long int: 33155719284104
    // lowres(/1) w original long long int: 0
    // 
    // Batu Ferrighi w original : 0


    /*
    
    // test code to examine the range for i ui li uli data type
    int n = 0;
    int i = n;
    unsigned int ui = n;
    long int li = n;
    unsigned long int uli = n;
    long long int lli = n;
    unsigned long long int ulli = n;
    cout << "sizeof(int)=" << sizeof(int) << ", sizeof(long int)=" << sizeof(long int) << ", sizeof (long long int)=" << sizeof(long long int) << endl;
    while (1) {
        i+=10;
        ui += 10;
        li += 10;
        uli += 10;
        lli += 10;
        ulli += 10;
        if (ulli % 100000000 == 0) {
            cout << i << '\t' << ui << '\t' << li << '\t' << uli << '\t' << lli<<'\t'<< ulli<<endl;
        }
    }
    // int and long int holding the same amount of range
    // unsigned int and unsigned long int holdings the same range

    *







    return 0;
}


*/


