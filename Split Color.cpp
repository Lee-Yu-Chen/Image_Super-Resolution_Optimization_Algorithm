
// modified from example in https://www.imagemagick.org/Magick++/Pixels.html

#include <Magick++.h> 
#include <iostream> 


using namespace std;
using namespace Magick;


int main(int argc, char** argv)
{
    InitializeMagick(*argv);

    Image image;

    // Read a file into image object 
    image.read("C:/Users/ACER/Downloads/Picture/Batu Ferringhi.jpg");

    string image_dir = "C:/Users/ACER/Desktop/folder/";     // folder needs to be create manually
    string filename_prefix = image_dir + ""; 
    // add some string to the filename to prevent overwritting the previous images

    
    image.write(filename_prefix+"original.png");



    // split image into RGB 3 images
    size_t width = image.size().width();
    size_t height = image.size().height();


	// create empty images
    Image red(Geometry(width, height), "#000000");
    Image green(Geometry(width, height), "#000000");
    Image blue(Geometry(width, height), "#000000");
    Image gray(Geometry(width, height), "#000000");
    Image cyan(Geometry(width, height), "#000000");
    Image magenta(Geometry(width, height), "#000000");
    Image yellow(Geometry(width, height), "#000000");

    red.type(TrueColorType);
    green.type(TrueColorType);
    blue.type(TrueColorType);
    gray.type(TrueColorType);
    cyan.type(TrueColorType);
    magenta.type(TrueColorType);
    yellow.type(TrueColorType);


    red.modifyImage();
    green.modifyImage();
    blue.modifyImage();
    gray.modifyImage();
    cyan.modifyImage();
    magenta.modifyImage();
    yellow.modifyImage();

    Pixels read(image);        // source image to read from
    Quantum* read_pix = read.get(0, 0, width, height);


    Pixels red_write(red);          // red image to generate
    Pixels green_write(green);          // green image to generate
    Pixels blue_write(blue);          // blue image to generate
    Pixels gray_write(gray);
    Pixels cyan_write(cyan);
    Pixels magenta_write(magenta);
    Pixels yellow_write(yellow);


    Quantum* red_write_p = red_write.get(0, 0, width, height);
    Quantum* green_write_p = green_write.get(0, 0, width, height);
    Quantum* blue_write_p = blue_write.get(0, 0, width, height);
    Quantum* gray_write_p = gray_write.get(0, 0, width, height);
    Quantum* cyan_write_p = cyan_write.get(0, 0, width, height);
    Quantum* magenta_write_p = magenta_write.get(0, 0, width, height);
    Quantum* yellow_write_p = yellow_write.get(0, 0, width, height);






	/*
	The data structure of the Quantum pointer
	the address pointed by the Quantum pointer is the red pixel value of the 1st pixel in the specified region
	the address pointed by the Quantum pointer + 1 is the green pixel value of the 1st pixel
	the address pointed by the Quantum pointer + 2 is the blue pixel value of the 1st pixel
	the address pointed by the Quantum pointer + 3 is the red pixel value of the 2nd pixel
	the address pointed by the Quantum pointer + 4 is the green pixel value of the 2nd pixel
	the address pointed by the Quantum pointer + 5 is the blue pixel value of the 2nd pixel
	... and so on
	
	R1 | G1 | B1 | R2 | G2 | B2 | ....
	
	
	*/



    int count = 0;

    for (int i = 0; i <= width - 1; i++) {

        for (int j = 0; j <= height - 1; j++) {
			
			// red image
            *(red_write_p + count) = *(read_pix + count);       // quantumRed()
            *(red_write_p + count + 1) = 0;                 // quantumGreen()
            *(red_write_p + count + 2) = 0;                 // quantumBlue()


			// green image
            *(green_write_p + count) = 0;
            *(green_write_p + count + 1) = *(read_pix + count + 1);
            *(green_write_p + count + 2) = 0;


			// blue image
            *(blue_write_p + count) = 0;
            *(blue_write_p + count + 1) = 0;
            *(blue_write_p + count + 2) = *(read_pix + count + 2);


			// gray image
            *(gray_write_p + count) = *(read_pix + count + 1);
            *(gray_write_p + count + 1) = *(read_pix + count + 1);
            *(gray_write_p + count + 2) = *(read_pix + count + 1);


			// cyan image
            *(cyan_write_p + count) = 0;
            *(cyan_write_p + count + 1) = *(read_pix + count + 1);
            *(cyan_write_p + count + 2) = *(read_pix + count + 1);


			// magenta image
            *(magenta_write_p + count) = *(read_pix + count + 1);
            *(magenta_write_p + count + 1) = 0;
            *(magenta_write_p + count + 2) = *(read_pix + count + 1);


			// yellow image
            *(yellow_write_p + count) = *(read_pix + count + 1);
            *(yellow_write_p + count + 1) = *(read_pix + count + 1);
            *(yellow_write_p + count + 2) = 0;



            count += 3;
        }
    }

    red_write.sync();
    red.write(filename_prefix+"red.png");

    green_write.sync();
    green.write(filename_prefix+"green.png");

    blue_write.sync();
    blue.write(filename_prefix+"blue.png");


    gray_write.sync();
    gray.write(filename_prefix+"gray.png");


    cyan_write.sync();
    cyan.write(filename_prefix+"cyan.png");


    magenta_write.sync();
    magenta.write(filename_prefix+"magenta.png");


    yellow_write.sync();
    yellow.write(filename_prefix+"yellow.png");

    return 0;
}


