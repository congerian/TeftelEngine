#pragma once 

#include <vector>
#include <functional>
#include <iostream>
#include <chrono>
#include <cmath>

#include <CL/cl.hpp>

#include <mingw_threading/mingw.shared_mutex.h>
#include <mingw_threading/mingw.thread.h>

#include "../public/particle.hxx"
#include "../public/util.hxx"

namespace TFN{class Engine;}

namespace TFN{

    class Engine{
        public:

        enum State{
            PAUSED = 0,
            PAUSING = 1,
            UNPAUSING = 2,
            WORKING = 3
        };

        enum Wall{
            LEFT = 0,
            TOP = 1, 
            RIGHT = 2,
            BOTTOM = 3
        };

        class IterationResult{
            public:
            
            Engine& engine; 


            FP dP;
            const unsigned long long deltaTime;

            IterationResult(Engine& _engine, FP _dP,  const unsigned long long& _deltaTime) : 
                dP{_dP},
                deltaTime{_deltaTime},
                engine{_engine}
            {}
        };

        class PhysicsEngine{
            public:

            Engine& engine;

            class CollisionResult{
                public: 
                
                using ParticleWallPair = std::pair < Particle*, Wall >;
                using ParticlePair = std::pair < Particle*, Particle* >;

                std::chrono::nanoseconds timeToCollision; 

                std::vector<ParticlePair> particleCollisions;
                std::vector<ParticleWallPair> wallCollisions;

                CollisionResult(std::chrono::nanoseconds _timeToCollision) : 
                    timeToCollision{_timeToCollision}
                {}

                CollisionResult(
                        std::chrono::nanoseconds _timeToCollision,
                        std::vector<ParticlePair> _particleCollisions,
                        std::vector<ParticleWallPair> _wallCollisions) :
                    timeToCollision{_timeToCollision},
                    particleCollisions{_particleCollisions},
                    wallCollisions{_wallCollisions}
                {}
            };



            CollisionResult calcClosestCollisions(std::chrono::nanoseconds deltaTime){

                std::chrono::nanoseconds minTime = deltaTime;

                std::vector<CollisionResult::ParticlePair> particleCollisions;
                std::vector<CollisionResult::ParticleWallPair> wallCollisions;

                //with other points

                Particle* it1;
                Particle* it2;

                size_t count = engine.nextParticles.size(); 

                it1 = engine.nextParticles.data();

                for(size_t i = 0; i < count; ++i){

                    FP vx = it1->vx;
                    FP vy = it1->vy;
                    FP x = it1->x;
                    FP y = it1->y;

                    FP R1 = it1->effectiveRadius;

                    if(vy > 1E-12){
                        std::chrono::nanoseconds t = std::chrono::nanoseconds((long long)(1E9 * (engine.top_wall_y - R1 - y) / vy));
                        if(t < minTime){
                            minTime = t;
                            wallCollisions.clear();
                            particleCollisions.clear();
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, TOP));
                        }
                        else if(t == minTime){
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, TOP));
                        }
                    }
                    else if(vy < -1E-12){
                        std::chrono::nanoseconds t = std::chrono::nanoseconds((long long)(1E9 * (engine.bottom_wall_y + R1 - y) / vy));
                        if(t < minTime){
                            minTime = t;
                            wallCollisions.clear();
                            particleCollisions.clear();
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, BOTTOM));
                        }
                        else if(t == minTime){
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, BOTTOM));
                        }
                    }
                    if(vx > 1E-12){
                        std::chrono::nanoseconds t = std::chrono::nanoseconds((long long)(1E9 * (engine.right_wall_x - R1 - x) / vx));
                        if(t < minTime){
                            minTime = t;
                            wallCollisions.clear();
                            particleCollisions.clear();
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, RIGHT));
                        }
                        else if(t == minTime){
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, RIGHT));
                        }
                    }
                    else if(vx < -1E-12){
                        std::chrono::nanoseconds t = std::chrono::nanoseconds((long long)(1E9 * (engine.left_wall_x + R1 - x) / vx));
                        if(t < minTime){
                            minTime = t;
                            wallCollisions.clear();
                            particleCollisions.clear();
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, LEFT));
                        }
                        else if(t == minTime){
                            wallCollisions.push_back(CollisionResult::ParticleWallPair(it1, LEFT));
                        }
                    }

                    it2 = it1;
                    ++it2;

                    for(size_t j = i + 1; j < count; ++j){
                        Particle& particleA = *it1;
                        Particle& particleB = *it2;

                        FP A_x = it1->x;
                        FP A_y = it1->y;

                        FP B_x = it1->vx - it2->vx;
                        FP B_y = it1->vy - it2->vy;

                        FP C_x = it2->x;
                        FP C_y = it2->y; 

                        FP r = it1->effectiveRadius + it2->effectiveRadius;

                        FP a = B_x*B_x + B_y*B_y;
                        
                        FP b = 2*(B_x*(A_x-C_x) + B_y*(A_y-C_y));

                        FP d = b*b - 4*a*((A_x-C_x)*(A_x-C_x) + (A_y-C_y)*(A_y-C_y) - r*r);

                        if (d < 0) {++it2; continue;};

                        FP e = - b - std::sqrt(d);

                        if(abs(a) < 1E-12) {++it2; continue;}

                        if (e > 0){
                            std::chrono::nanoseconds t = std::chrono::nanoseconds((long long)(1E9 * e / (2*a)));
                            if (t < minTime){
                                particleCollisions.clear();
                                wallCollisions.clear();
                                particleCollisions.push_back(CollisionResult::ParticlePair(it1, it2));
                                minTime = t;
                            }
                            ++it2;
                            continue;
                        } 

                        e = -2*b - 3;

                        if (e > 0){
                            std::chrono::nanoseconds t = std::chrono::nanoseconds((long long)(1E9 * e / (2*a)));
                            if (t < minTime){
                                particleCollisions.clear();
                                wallCollisions.clear();
                                particleCollisions.push_back(CollisionResult::ParticlePair(it1, it2));
                                minTime = t;
                            }
                        } 

                        ++it2;
                    }

                    ++it1;
                }


                //with walls


                return CollisionResult(minTime, particleCollisions, wallCollisions);
            }

            void linearMove(std::chrono::nanoseconds deltaTime){

                Particle* it;

                size_t count = engine.nextParticles.size(); 

                it = engine.nextParticles.data();

                for(size_t i = 0; i < count; ++i){
                    it->x = it->x + it->vx * (deltaTime.count() * 1E-9);
                    it->y = it->y + it->vy * (deltaTime.count() * 1E-9);

                    ++it;
                }

                return;
            }

            FP processCollisions(CollisionResult collisionResult){
                
                FP dP = 0.0;

                CollisionResult::ParticlePair* it1;
                CollisionResult::ParticleWallPair* it2;

                size_t particlePairCount = collisionResult.particleCollisions.size(); 
                size_t particleWallPairCount = collisionResult.wallCollisions.size(); 

                it1 = collisionResult.particleCollisions.data();
                it2 = collisionResult.wallCollisions.data();
                
                for(size_t i = 0; i < particlePairCount; ++i){
                    Particle* p1 = it1->first;
                    Particle* p2 = it1->second;

                    FP r1_x = p1->x;
                    FP r1_y = p1->y;

                    FP r2_x = p2->x;
                    FP r2_y = p2->y;
                    
                    FP v1_x = p1->vx;
                    FP v1_y = p1->vy;
                    
                    FP v2_x = p2->vx;
                    FP v2_y = p2->vy;
                    
                    FP R1 = p1->effectiveRadius;
                    FP R2 = p2->effectiveRadius;
                    
                    FP m1 = p1->mass;
                    FP m2 = p2->mass;



                    FP e_a_x = (r2_x - r1_x)/(R1+R2);
                    FP e_a_y = (r2_y - r1_y)/(R1+R2);

                    FP e_b_x = e_a_y;
                    FP e_b_y = -e_a_x;

                    FP v1_a = v1_x*e_a_x + v1_y*e_a_y;
                    FP v2_a = v2_x*e_a_x + v2_y*e_a_y;

                    FP v1_b = v1_x*e_b_x + v1_y*e_b_y;
                    FP v2_b = v2_x*e_b_x + v2_y*e_b_y;

                    FP new_v1_a = (2*m2*v2_a + v1_a*(m1-m2))/(m1 + m2);
                    FP new_v2_a = (2*m1*v1_a + v2_a*(m2-m1))/(m1 + m2);



                    p1->vx = new_v1_a * e_a_x + v1_b * e_b_x;
                    p1->vy = new_v1_a * e_a_y + v1_b * e_b_y;

                    p2->vx = new_v2_a * e_a_x + v2_b * e_b_x;
                    p2->vy = new_v2_a * e_a_y + v2_b * e_b_y;

                    ++it1;
                }

                for(size_t i = 0; i < particleWallPairCount; ++i){
                    Particle* p = it2->first;
                    
                    FP vx = p->vx;
                    FP vy = p->vy;

                    Wall w = it2->second;
                    FP m = p->mass;

                    if(w == TOP && vy > 1E-12) {vy = -vy; dP -= m*2*vy;} 
                    if(w == BOTTOM && vy < -1E-12) {vy = -vy; dP += m*2*vy;}
                    if(w == RIGHT && vx > 1E-12) {vx = -vx; dP -= m*2*vx;}
                    if(w == LEFT && vx < -1E-12) {vx = -vx; dP += m*2*vx;}

                    p->vx = vx;
                    p->vy = vy;

                    ++it2;
                }

                return dP;
            }

            FP processPhysics(std::chrono::nanoseconds deltaTime){
                
                std::shared_lock shared_lock(engine.nextParticlesMutex, std::defer_lock);
                std::unique_lock unique_lock(engine.nextParticlesMutex, std::defer_lock);

                FP dP = 0.0;

                if (deltaTime.count() > 0){
                    using namespace std::chrono_literals;
                    

                    shared_lock.lock();
                    CollisionResult collisionResult = calcClosestCollisions(deltaTime);
                    shared_lock.unlock();

                    unique_lock.lock();
                    linearMove(collisionResult.timeToCollision);
                    unique_lock.unlock();

                    unique_lock.lock();
                    dP += processCollisions(collisionResult);
                    unique_lock.unlock();      

                    processPhysics(deltaTime - collisionResult.timeToCollision); 
                }

                return dP;
            }

            PhysicsEngine(Engine& _engine) :
                engine{_engine}
            {}
        };

        //Getters:

        State getState();
        std::vector<Particle> getCurrentParticleState()
        {
            std::shared_lock shared_lock(currentParticlesMutex, std::defer_lock);

            shared_lock.lock();
            std::vector<Particle> copiedParticles = currentParticles;
            shared_lock.unlock();

            return copiedParticles;
        }

        //Modifying particles:

        void deleteAllParticles();
        void spawnParticle(Particle);
        void spawnParticles(std::vector<Particle> request){
            std::unique_lock unique_lock(queueMutex, std::defer_lock);
            
            unique_lock.lock();
            spawnQueue.insert(spawnQueue.end(), request.begin(), request.end());
            unique_lock.unlock();
        }

        //Callback

        std::function<void(IterationResult)> onIter;

        //Constructor

        Engine() : 
            physicsEngine{PhysicsEngine(*this)}    
        {
            mainThreadId = std::this_thread::get_id();
            workingThread = std::thread(&Engine::workingLoop, this);
            workingThreadId = workingThread.get_id();

        }

        std::thread::id mainThreadId;
        std::thread::id workingThreadId;

        FP top_wall_y = 10;
        FP right_wall_x = 10;
        FP bottom_wall_y = -10;
        FP left_wall_x = -10;

        private:

        std::shared_mutex currentParticlesMutex;
        std::shared_mutex nextParticlesMutex;
        std::shared_mutex queueMutex;

        PhysicsEngine physicsEngine;

        std::vector<Particle> spawnQueue;

        void workingLoop(){
            using namespace std::chrono_literals;

            std::cout << "Starting working loop!" << std::endl;

            while(1){
                if(state == PAUSING) state = PAUSED;
                if(state == UNPAUSING) state = WORKING;

                //Spawning new particles

                std::shared_lock nextParticlesSharedLock(nextParticlesMutex, std::defer_lock);
                std::unique_lock nextParticlesUniqueLock(nextParticlesMutex, std::defer_lock);
                std::shared_lock currentParticlesSharedLock(currentParticlesMutex, std::defer_lock);
                std::unique_lock currentParticlesUniqueLock(currentParticlesMutex, std::defer_lock);
                std::shared_lock queueSharedLock(queueMutex, std::defer_lock);
                std::unique_lock queueUniqueLock(queueMutex, std::defer_lock);


                std::lock(nextParticlesUniqueLock, queueUniqueLock);                

                nextParticles.insert(nextParticles.end(), spawnQueue.begin(), spawnQueue.end());
                
                nextParticlesUniqueLock.unlock();

                spawnQueue.clear();
                
                queueUniqueLock.unlock();


                //physics      

                auto start = std::chrono::high_resolution_clock::now();
                
                FP dP = physicsEngine.processPhysics(std::chrono::nanoseconds(250));


                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<FP, std::milli> elapsed = end-start;

                //graphics

                

                //callback


                std::lock(currentParticlesUniqueLock, nextParticlesSharedLock);

                currentParticles = nextParticles;

                nextParticlesSharedLock.unlock();
                currentParticlesUniqueLock.unlock();

                //currentParticlesUniqueLock.lock();

                if(onIter){
                    IterationResult iterationResult(*this, dP, elapsed.count());
                    onIter(iterationResult);
                }

                //currentParticlesUniqueLock.unlock();

                

            }

        }

        std::thread workingThread;

        std::vector<Particle> currentParticles;
        std::vector<Particle> nextParticles;

        State state;

    };
} 

