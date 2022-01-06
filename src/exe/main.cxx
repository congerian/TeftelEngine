#include "main.hxx"

#include "particle.hxx"
#include "../public/util.hxx"


#include <random>
#include <iomanip>
#include <fstream>
#include <iostream>

#include <SDL2/SDL.h>
#include <CL/cl.hpp>

int HEIGHT = 600;
int WIDTH = 800;

TFN::FP scaling = 1E-1;

TFN::FP widthMeters = WIDTH * scaling;
TFN::FP heightMeters = HEIGHT * scaling;

unsigned long long my_time = 0;
TFN::FP P = 0.0;

void DrawCircle(SDL_Renderer * renderer, int32_t centreX, int32_t centreY, int32_t radius)
{
   const int32_t diameter = (radius * 2);

   int32_t x = (radius - 1);
   int32_t y = 0;
   int32_t tx = 1;
   int32_t ty = 1;
   int32_t error = (tx - diameter);

   while (x >= y)
   {
      SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
      SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
      SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
      SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
      SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
      SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
      SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
      SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);

      if (error <= 0)
      {
         ++y;
         error += ty;
         ty += 2;
      }

      if (error > 0)
      {
         --x;
         tx += 2;
         error += (tx - diameter);
      }
   }
}

void DrawParticle(SDL_Renderer * renderer, TFN::Particle particle){
	int x = WIDTH/2 + particle.x / scaling;
	int y = HEIGHT/2 + particle.y / scaling;
	int r = particle.effectiveRadius / scaling;

	DrawCircle(renderer, x, y, r);
}

int main()
{
	TFN::Engine engine;

	engine.left_wall_x = -widthMeters/2;
	engine.right_wall_x = widthMeters/2;
	engine.top_wall_y = heightMeters/2;
	engine.bottom_wall_y = -heightMeters/2;

	std::vector<TFN::Particle> newParticles;

	std::random_device rd{};
    std::mt19937 gen{rd()};

	std::uniform_real_distribution<> genX(-widthMeters/6, widthMeters/6);
	std::uniform_real_distribution<> genY(-heightMeters/3, heightMeters/3);
	
	std::uniform_real_distribution<> genV(-1E1, 1E1);

	/*
	for(int i = 0; i < 200; ++i){
		newParticles.push_back(TFN::Particle(genX(gen), genY(gen), genV(gen), genV(gen), 3E-10, 7E-1));
	}
	*/

	/*
	for(int i = 0; i < 300; ++i){
		newParticles.push_back(TFN::Particle(genX(gen), genY(gen), 1E6, 1E6, 1E-10, 3E-1));
	}
	*/

	for(int i = -8; i < 8; ++i){
		for(int j = -15; j < 15; ++j){
			if(j % 2 == 0)
				newParticles.push_back(
					TFN::Particle(
						(widthMeters/2/20)*i, (heightMeters/2/20)*j, genV(gen), genV(gen), 3E-10, 7E-1));
			else
				newParticles.push_back(
					TFN::Particle(
						(widthMeters/2/20)*i, (heightMeters/2/20)*j, genV(gen), genV(gen), 1E-10, 3E-1));
		}
	}


    
	newParticles.push_back(TFN::Particle(-widthMeters/3, 0, -1E6, 0, 1E-7, 1E1));

	engine.spawnParticles(newParticles);

	std::ofstream dataFile;
	dataFile.open("data.csv", std::ios::trunc);
	dataFile.close();

	int iterationNumber = 0;

	std::function<void(TFN::Engine::IterationResult)> callback = [&](TFN::Engine::IterationResult iterationResult){
				++iterationNumber;

				dataFile.open("data.csv", std::ios::out | std::ios::app);

				P = (P*my_time + iterationResult.dP * iterationResult.deltaTime)/(my_time+iterationResult.deltaTime);
				my_time += iterationResult.deltaTime;

                std::cout 	<< "Iteration took: " << std::setw(3) << iterationResult.deltaTime 
							<< "ms! dP = " << std::setw(12) << iterationResult.dP << " N/M.   " 
							<< "P_avg = " << std::setw(12) << P << " N/M. \r";
				std::cout.flush();

				size_t count = engine.getCurrentParticleState().size(); 

                TFN::Particle* it = engine.getCurrentParticleState().data();

				std::stringstream str;

				str << iterationNumber;

                for(size_t i = 0; i < count; ++i){
					TFN::Particle p = *it;
					str << "," << std::to_string(p.mass * ((p.vx)*(p.vx) + (p.vy)*(p.vy))/2);
					++it;
				}

				str << "\n";

				dataFile << str.str();
				dataFile.close();
				
				return;
            };

	engine.onIter = callback;

	SDL_Init(SDL_INIT_VIDEO);

	SDL_Window *window = SDL_CreateWindow(
		"SDL2Test",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		WIDTH,
		HEIGHT,
		SDL_WINDOW_VULKAN
	);

	SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
	SDL_Event e;


	while(1){
		while (SDL_PollEvent(&e) != 0)
		{
			if(e.type == SDL_QUIT) return 0;		
		}

		SDL_SetRenderDrawColor(renderer, 0xFF, 0xFF, 0xFF, 0xFF);
		SDL_RenderClear(renderer);
		
		SDL_SetRenderDrawColor(renderer, 255, 0, 0, 0);
		
		std::vector<TFN::Particle> particles = engine.getCurrentParticleState();
		for(auto p : particles){
			DrawParticle(renderer, p);
		}
		SDL_RenderPresent(renderer);

	}

	while(1){
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(1000ms);
		std::cout << "I am aLAIVe!" << std::endl;
	}
	return 0;
}

