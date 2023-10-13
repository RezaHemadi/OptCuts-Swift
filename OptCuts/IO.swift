//
//  IO.swift
//  OptCuts
//
//  Created by Reza on 10/13/23.
//

import Foundation
import ModelIO
import MetalKit
import Matrix

/// Travese scene graph rooted at object and run a closure on each element,
///  passing on each element's flattened node index as well as it's parent index
func walkSceneGraph(rootAt object: MDLObject, perNodeBody: (MDLObject, Int, Int?) -> Void) {
    var currentIndex = 0
    
    func walkGraph(object: MDLObject, currentIndex: inout Int, parentIndex: Int?, perNodeBody: (MDLObject, Int, Int?) -> Void) {
        perNodeBody(object, currentIndex, parentIndex)
        
        let ourIndex = currentIndex
        currentIndex += 1
        var pointer: ObjCBool = false
        var count: Int = 0
        object.enumerateChildObjects(of: MDLObject.self, root: object, using: { _, _ in count += 1 }, stopPointer: &pointer)
        guard count != 0 else { return }
        
        for childIndex in 0..<object.children.count {
            walkGraph(
                object: object.children[childIndex],
                currentIndex: &currentIndex,
                parentIndex: ourIndex,
                perNodeBody: perNodeBody
            )
        }
    }
    
    walkGraph(object: object, currentIndex: &currentIndex, parentIndex: nil, perNodeBody: perNodeBody)
}


class GeometryLoader {
    static func LoadMesh(allocator: MTKMeshBufferAllocator,
                         device: MTLDevice,
                         url: URL,
                         vertices: inout Mat<Float>,
                         normals: inout Mat<Float>,
                         faces: inout Mat<UInt32>) {
        let mtlVertexDescriptor = MTLVertexDescriptor()
        mtlVertexDescriptor.attributes[0].format = .float3
        mtlVertexDescriptor.attributes[0].offset = 0
        mtlVertexDescriptor.attributes[0].bufferIndex = 0
        mtlVertexDescriptor.attributes[1].format = .float3
        mtlVertexDescriptor.attributes[1].offset = 3 * MemoryLayout<Float>.size
        mtlVertexDescriptor.attributes[1].bufferIndex = 0
        mtlVertexDescriptor.layouts[0].stride = MemoryLayout<Float>.stride * 6
        
        let vertexDesc = MTKModelIOVertexDescriptorFromMetal(mtlVertexDescriptor)
        (vertexDesc.attributes[0] as! MDLVertexAttribute).name = MDLVertexAttributePosition
        (vertexDesc.attributes[1] as! MDLVertexAttribute).name = MDLVertexAttributeNormal
        
        let asset = MDLAsset(url: url, vertexDescriptor: vertexDesc, bufferAllocator: allocator)
        let root = asset.object(at: 0)
        
        var localTransforms: [simd_float4x4] = []
        var meshTransforms: [simd_float4x4] = []
        var parentIndices: [Int?] = []
        
        var meshNodeIndices: [Int] = []
        var meshPositionsBuffers: [MTLBuffer] = []
        var meshPositionsBuffersOffsets: [Int] = []
        var meshVertexCounts: [Int] = []
        var indexCount: [[Int]] = [[]]
        var indexType: [[MTLIndexType]] = [[]]
        var indexBuffers: [[MTLBuffer]] = [[]]
        var indexBuffersOffsets: [[Int]] = [[]]
        var primitiveTypes: [[MTLPrimitiveType]] = [[]]
        
        walkSceneGraph(rootAt: root) { object, curIdx, _ in
            if let mesh = object as? MDLMesh {
                assert(mesh.vertexBuffers.count == 1)
                
                do {
                    meshNodeIndices.append(curIdx)
                    let mtkMesh = try MTKMesh(mesh: mesh, device: device)
                    meshPositionsBuffers.append(mtkMesh.vertexBuffers[0].buffer)
                    
                    
                    meshPositionsBuffersOffsets.append(mtkMesh.vertexBuffers[0].offset)
                    meshVertexCounts.append(mesh.vertexCount)
                    
                    indexCount.append([])
                    indexType.append([])
                    indexBuffers.append([])
                    indexBuffersOffsets.append([])
                    primitiveTypes.append([])
                    
                    for i in 0..<mtkMesh.submeshes.count {
                        let subMesh = mtkMesh.submeshes[i]
                        indexCount[meshPositionsBuffers.count - 1].append(subMesh.indexCount)
                        indexType[meshPositionsBuffers.count - 1].append(subMesh.indexType)
                        indexBuffers[meshPositionsBuffers.count - 1].append(subMesh.indexBuffer.buffer)
                        indexBuffersOffsets[meshPositionsBuffers.count - 1].append(subMesh.indexBuffer.offset)
                        primitiveTypes[meshPositionsBuffers.count - 1].append(subMesh.primitiveType)
                    }
                } catch {
                    print("error creating mesh \(error)")
                }
            }
        }
        
        // flatten hierarchy
        walkSceneGraph(rootAt: root) { object, currentIdx, parentIndx in
            if let transform = object.transform {
                localTransforms.append(transform.matrix)
            } else {
                localTransforms.append(matrix_identity_float4x4)
            }
            parentIndices.append(parentIndx)
        }
        
        // calcualte the transform of each mesh relative to root and store in meshtransforms
        for meshNodeIndex in 0..<meshNodeIndices.count {
            let index = meshNodeIndices[meshNodeIndex]
            var parent = parentIndices[index]
            var myTransform = localTransforms[index]
            
            while parent != nil {
                myTransform = localTransforms[parent!] * myTransform
                parent = parentIndices[parent!]
            }
            
            meshTransforms.append(myTransform)
        }
        
        // store vertices
        let totalVertexCount: Int = meshVertexCounts.reduce(0, { $0 + $1})
        // resize input vertices matrix
        vertices.resize(totalVertexCount, 3)
        normals.resize(totalVertexCount, 3)
        vertices.setZero()
        var processedVertices: Int = 0
        for meshNodeIndex in 0..<meshNodeIndices.count {
            let modelSpaceTransform = meshTransforms[meshNodeIndex]
            
            for i in 0..<meshVertexCounts[meshNodeIndex] {
                let xOffset = MemoryLayout<Float>.size * (6 * i)
                let yOffset = MemoryLayout<Float>.size * (6 * i + 1)
                let zOffset = MemoryLayout<Float>.size * (6 * i + 2)
                
                let normalXOffset = MemoryLayout<Float>.size * (6 * i + 3)
                let normalYOffset = MemoryLayout<Float>.size * (6 * i + 4)
                let normalZOffset = MemoryLayout<Float>.size * (6 * i + 5)
                
                let x: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + xOffset).assumingMemoryBound(to: Float.self).pointee
                let y: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + yOffset).assumingMemoryBound(to: Float.self).pointee
                let z: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + zOffset).assumingMemoryBound(to: Float.self).pointee
                
                let normalX: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + normalXOffset).assumingMemoryBound(to: Float.self).pointee
                let normalY: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + normalYOffset).assumingMemoryBound(to: Float.self).pointee
                let normalZ: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + normalZOffset).assumingMemoryBound(to: Float.self).pointee
                
                // get the point in model space
                let localSpaceHPoint: SIMD4<Float> = .init(x: x, y: y, z: z, w: 1.0)
                let modelSpaceHPoint = modelSpaceTransform * localSpaceHPoint
                
                // get the normal in model space
                let localSpaceHNormal: SIMD4<Float> = .init(x: normalX, y: normalY, z: normalZ, w: 0.0)
                let modelSpaceHNormal = modelSpaceTransform * localSpaceHNormal
                
                vertices[processedVertices, 0] = modelSpaceHPoint.x
                vertices[processedVertices, 1] = modelSpaceHPoint.y
                vertices[processedVertices, 2] = modelSpaceHPoint.z
                
                normals[processedVertices, 0] = modelSpaceHNormal.x
                normals[processedVertices, 1] = modelSpaceHNormal.y
                normals[processedVertices, 2] = modelSpaceHNormal.z
                
                processedVertices += 1
            }
            
            let meshIndexBuffers = indexBuffers[meshNodeIndex]
            let meshIndexCounts = indexCount[meshNodeIndex]
            
            for i in 0..<meshIndexBuffers.count {
                let indexBuffer = meshIndexBuffers[i]
                let indexCount = meshIndexCounts[i]
                
                assert(indexCount % 3 == 0)
                let triangleCount = indexCount / 3
                faces.resize(triangleCount, 3)
                let indexBufferAddress = indexBuffer.contents().advanced(by: indexBuffersOffsets[meshNodeIndex][i])
                
                // loop over triangles
                for t in 0..<triangleCount {
                    (faces.valuesPtr.pointer + (3 * t)).initialize(
                        from: indexBufferAddress.assumingMemoryBound(to: UInt32.self).advanced(by: 3 * t),
                        count: 3)
                }
            }
        }
    }
    
    static func LoadMesh(allocator: MTKMeshBufferAllocator,
                         device: MTLDevice,
                         url: URL,
                         vertices: inout Mat<Float>,
                         faces: inout Mat<UInt32>) {
        let mtlVertexDescriptor = MTLVertexDescriptor()
        mtlVertexDescriptor.attributes[0].format = .float3
        mtlVertexDescriptor.attributes[0].offset = 0
        mtlVertexDescriptor.attributes[0].bufferIndex = 0
        mtlVertexDescriptor.layouts[0].stride = MemoryLayout<Float>.stride * 3
        
        let vertexDesc = MTKModelIOVertexDescriptorFromMetal(mtlVertexDescriptor)
        (vertexDesc.attributes[0] as! MDLVertexAttribute).name = MDLVertexAttributePosition
        
        let asset = MDLAsset(url: url, vertexDescriptor: vertexDesc, bufferAllocator: allocator)
        let root = asset.object(at: 0)
        
        var localTransforms: [simd_float4x4] = []
        var meshTransforms: [simd_float4x4] = []
        var parentIndices: [Int?] = []
        
        var meshNodeIndices: [Int] = []
        var meshPositionsBuffers: [MTLBuffer] = []
        var meshPositionsBuffersOffsets: [Int] = []
        var meshVertexCounts: [Int] = []
        var indexCount: [[Int]] = [[]]
        var indexType: [[MTLIndexType]] = [[]]
        var indexBuffers: [[MTLBuffer]] = [[]]
        var indexBuffersOffsets: [[Int]] = [[]]
        var primitiveTypes: [[MTLPrimitiveType]] = [[]]
        
        walkSceneGraph(rootAt: root) { object, curIdx, _ in
            if let mesh = object as? MDLMesh {
                assert(mesh.vertexBuffers.count == 1)
                
                do {
                    meshNodeIndices.append(curIdx)
                    let mtkMesh = try MTKMesh(mesh: mesh, device: device)
                    meshPositionsBuffers.append(mtkMesh.vertexBuffers[0].buffer)
                    
                    
                    meshPositionsBuffersOffsets.append(mtkMesh.vertexBuffers[0].offset)
                    meshVertexCounts.append(mesh.vertexCount)
                    
                    indexCount.append([])
                    indexType.append([])
                    indexBuffers.append([])
                    indexBuffersOffsets.append([])
                    primitiveTypes.append([])
                    
                    for i in 0..<mtkMesh.submeshes.count {
                        let subMesh = mtkMesh.submeshes[i]
                        indexCount[meshPositionsBuffers.count - 1].append(subMesh.indexCount)
                        indexType[meshPositionsBuffers.count - 1].append(subMesh.indexType)
                        indexBuffers[meshPositionsBuffers.count - 1].append(subMesh.indexBuffer.buffer)
                        indexBuffersOffsets[meshPositionsBuffers.count - 1].append(subMesh.indexBuffer.offset)
                        primitiveTypes[meshPositionsBuffers.count - 1].append(subMesh.primitiveType)
                    }
                } catch {
                    print("error creating mesh \(error)")
                }
            }
        }
        
        // flatten hierarchy
        walkSceneGraph(rootAt: root) { object, currentIdx, parentIndx in
            if let transform = object.transform {
                localTransforms.append(transform.matrix)
            } else {
                localTransforms.append(matrix_identity_float4x4)
            }
            parentIndices.append(parentIndx)
        }
        
        // calcualte the transform of each mesh relative to root and store in meshtransforms
        for meshNodeIndex in 0..<meshNodeIndices.count {
            let index = meshNodeIndices[meshNodeIndex]
            var parent = parentIndices[index]
            var myTransform = localTransforms[index]
            
            while parent != nil {
                myTransform = localTransforms[parent!] * myTransform
                parent = parentIndices[parent!]
            }
            
            meshTransforms.append(myTransform)
        }
        
        // store vertices
        let totalVertexCount: Int = meshVertexCounts.reduce(0, { $0 + $1})
        // resize input vertices matrix
        vertices.resize(totalVertexCount, 3)
        vertices.setZero()
        var processedVertices: Int = 0
        for meshNodeIndex in 0..<meshNodeIndices.count {
            let modelSpaceTransform = meshTransforms[meshNodeIndex]
            
            for i in 0..<meshVertexCounts[meshNodeIndex] {
                let xOffset = MemoryLayout<Float>.size * (3 * i)
                let yOffset = MemoryLayout<Float>.size * (3 * i + 1)
                let zOffset = MemoryLayout<Float>.size * (3 * i + 2)
                
                let x: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + xOffset).assumingMemoryBound(to: Float.self).pointee
                let y: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + yOffset).assumingMemoryBound(to: Float.self).pointee
                let z: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + zOffset).assumingMemoryBound(to: Float.self).pointee
                
                // get the point in model space
                let localSpaceHPoint: SIMD4<Float> = .init(x: x, y: y, z: z, w: 1.0)
                let modelSpaceHPoint = modelSpaceTransform * localSpaceHPoint
                
                vertices[processedVertices, 0] = modelSpaceHPoint.x
                vertices[processedVertices, 1] = modelSpaceHPoint.y
                vertices[processedVertices, 2] = modelSpaceHPoint.z
                processedVertices += 1
            }
            
            let meshIndexBuffers = indexBuffers[meshNodeIndex]
            let meshIndexCounts = indexCount[meshNodeIndex]
            
            for i in 0..<meshIndexBuffers.count {
                let indexBuffer = meshIndexBuffers[i]
                let indexCount = meshIndexCounts[i]
                
                assert(indexCount % 3 == 0)
                let triangleCount = indexCount / 3
                faces.resize(triangleCount, 3)
                let indexBufferAddress = indexBuffer.contents().advanced(by: indexBuffersOffsets[meshNodeIndex][i])
                
                // loop over triangles
                for t in 0..<triangleCount {
                    (faces.valuesPtr.pointer + (3 * t)).initialize(
                        from: indexBufferAddress.assumingMemoryBound(to: UInt32.self).advanced(by: 3 * t),
                        count: 3)
                }
            }
        }
    }
    
    static func WriteToFile(_ url: URL,
                            allocator: MTKMeshBufferAllocator,
                            vertexData: Mat<Float>,
                            faceData: Mat<UInt32>,
                            uvData: Mat<Float>? = nil,
                            normalData: Mat<Float>? = nil) throws {
        let asset = MDLAsset()
        let data = Data.init(bytes: vertexData.valuesPtr.pointer,
                             count: MemoryLayout<Float>.size * vertexData.size.count)
        let vBuffer = allocator.newBuffer(with: data, type: .vertex)
        
        let indexData = Data.init(bytes: faceData.valuesPtr.pointer,
                                  count: MemoryLayout<UInt32>.size * faceData.size.count)
        let iBuffer = allocator.newBuffer(with: indexData, type: .index)
        let submesh = MDLSubmesh(indexBuffer: iBuffer,
                                 indexCount: faceData.size.count,
                                 indexType: .uint32,
                                 geometryType: .triangles,
                                 material: nil)
        
        var vBuffers: [MDLMeshBuffer] = [vBuffer]
        
        let vDescriptor = MDLVertexDescriptor()
        vDescriptor.attributes[0] = MDLVertexAttribute(name: MDLVertexAttributePosition,
                                                       format: .float3,
                                                       offset: 0,
                                                       bufferIndex: 0)
        vDescriptor.layouts[0] = MDLVertexBufferLayout(stride: MemoryLayout<Float>.size * 3)
        
        
        if let uvData = uvData {
            let uvData = Data.init(bytes: uvData.valuesPtr.pointer,
                                   count: MemoryLayout<Float>.size * uvData.size.count)
            let uvBuffer = allocator.newBuffer(with: uvData, type: .vertex)
            vBuffers.append(uvBuffer)
            
            vDescriptor.attributes[vBuffers.count - 1] = MDLVertexAttribute(name: MDLVertexAttributeTextureCoordinate,
                                                           format: .float2,
                                                           offset: 0,
                                                           bufferIndex: vBuffers.count - 1)
            vDescriptor.layouts[vBuffers.count - 1] = MDLVertexBufferLayout(stride: MemoryLayout<Float>.size * 2)
        }
        
        if let normalData = normalData {
            let normalData = Data.init(bytes: normalData.valuesPtr.pointer,
                                       count: MemoryLayout<Float>.size * normalData.size.count)
            let normalBuffer = allocator.newBuffer(with: normalData, type: .vertex)
            vBuffers.append(normalBuffer)
            
            vDescriptor.attributes[vBuffers.count - 1] = MDLVertexAttribute(name: MDLVertexAttributeNormal,
                                                           format: .float3,
                                                           offset: 0,
                                                           bufferIndex: vBuffers.count - 1)
            vDescriptor.layouts[vBuffers.count - 1] = MDLVertexBufferLayout(stride: MemoryLayout<Float>.size * 3)
        }
        
        let mdlMesh = MDLMesh(vertexBuffers: vBuffers,
                              vertexCount: vertexData.rows,
                              descriptor: vDescriptor,
                              submeshes: [submesh])
        asset.add(mdlMesh)
        
        FileManager.default.createFile(atPath: url.path(), contents: nil)
        do {
            try asset.export(to: url)
        }
    }
    
    static func LoadMesh(allocator: MTKMeshBufferAllocator,
                         device: MTLDevice,
                         url: URL,
                         vertices: inout Mat<Double>,
                         faces: inout Mat<Int>) {
        let mtlVertexDescriptor = MTLVertexDescriptor()
        mtlVertexDescriptor.attributes[0].format = .float3
        mtlVertexDescriptor.attributes[0].offset = 0
        mtlVertexDescriptor.attributes[0].bufferIndex = 0
        mtlVertexDescriptor.layouts[0].stride = MemoryLayout<Float>.stride * 3
        
        let vertexDesc = MTKModelIOVertexDescriptorFromMetal(mtlVertexDescriptor)
        (vertexDesc.attributes[0] as! MDLVertexAttribute).name = MDLVertexAttributePosition
        
        let asset = MDLAsset(url: url, vertexDescriptor: vertexDesc, bufferAllocator: allocator)
        let root = asset.object(at: 0)
        
        var localTransforms: [simd_float4x4] = []
        var meshTransforms: [simd_float4x4] = []
        var parentIndices: [Int?] = []
        
        var meshNodeIndices: [Int] = []
        var meshPositionsBuffers: [MTLBuffer] = []
        var meshPositionsBuffersOffsets: [Int] = []
        var meshVertexCounts: [Int] = []
        var indexCount: [[Int]] = [[]]
        var indexType: [[MTLIndexType]] = [[]]
        var indexBuffers: [[MTLBuffer]] = [[]]
        var indexBuffersOffsets: [[Int]] = [[]]
        var primitiveTypes: [[MTLPrimitiveType]] = [[]]
        
        walkSceneGraph(rootAt: root) { object, curIdx, _ in
            if let mesh = object as? MDLMesh {
                assert(mesh.vertexBuffers.count == 1)
                
                do {
                    meshNodeIndices.append(curIdx)
                    let mtkMesh = try MTKMesh(mesh: mesh, device: device)
                    meshPositionsBuffers.append(mtkMesh.vertexBuffers[0].buffer)
                    
                    
                    meshPositionsBuffersOffsets.append(mtkMesh.vertexBuffers[0].offset)
                    meshVertexCounts.append(mesh.vertexCount)
                    
                    indexCount.append([])
                    indexType.append([])
                    indexBuffers.append([])
                    indexBuffersOffsets.append([])
                    primitiveTypes.append([])
                    
                    for i in 0..<mtkMesh.submeshes.count {
                        let subMesh = mtkMesh.submeshes[i]
                        indexCount[meshPositionsBuffers.count - 1].append(subMesh.indexCount)
                        indexType[meshPositionsBuffers.count - 1].append(subMesh.indexType)
                        indexBuffers[meshPositionsBuffers.count - 1].append(subMesh.indexBuffer.buffer)
                        indexBuffersOffsets[meshPositionsBuffers.count - 1].append(subMesh.indexBuffer.offset)
                        primitiveTypes[meshPositionsBuffers.count - 1].append(subMesh.primitiveType)
                    }
                } catch {
                    print("error creating mesh \(error)")
                }
            }
        }
        
        // flatten hierarchy
        walkSceneGraph(rootAt: root) { object, currentIdx, parentIndx in
            if let transform = object.transform {
                localTransforms.append(transform.matrix)
            } else {
                localTransforms.append(matrix_identity_float4x4)
            }
            parentIndices.append(parentIndx)
        }
        
        // calcualte the transform of each mesh relative to root and store in meshtransforms
        for meshNodeIndex in 0..<meshNodeIndices.count {
            let index = meshNodeIndices[meshNodeIndex]
            var parent = parentIndices[index]
            var myTransform = localTransforms[index]
            
            while parent != nil {
                myTransform = localTransforms[parent!] * myTransform
                parent = parentIndices[parent!]
            }
            
            meshTransforms.append(myTransform)
        }
        
        // store vertices
        let totalVertexCount: Int = meshVertexCounts.reduce(0, { $0 + $1})
        // resize input vertices matrix
        vertices.resize(totalVertexCount, 3)
        vertices.setZero()
        var processedVertices: Int = 0
        for meshNodeIndex in 0..<meshNodeIndices.count {
            let modelSpaceTransform = meshTransforms[meshNodeIndex]
            
            for i in 0..<meshVertexCounts[meshNodeIndex] {
                let xOffset = MemoryLayout<Float>.size * (3 * i)
                let yOffset = MemoryLayout<Float>.size * (3 * i + 1)
                let zOffset = MemoryLayout<Float>.size * (3 * i + 2)
                
                let x: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + xOffset).assumingMemoryBound(to: Float.self).pointee
                let y: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + yOffset).assumingMemoryBound(to: Float.self).pointee
                let z: Float = meshPositionsBuffers[meshNodeIndex].contents().advanced(by: meshPositionsBuffersOffsets[meshNodeIndex] + zOffset).assumingMemoryBound(to: Float.self).pointee
                
                // get the point in model space
                let localSpaceHPoint: SIMD4<Float> = .init(x: x, y: y, z: z, w: 1.0)
                let modelSpaceHPoint = modelSpaceTransform * localSpaceHPoint
                
                vertices[processedVertices, 0] = Double(modelSpaceHPoint.x)
                vertices[processedVertices, 1] = Double(modelSpaceHPoint.y)
                vertices[processedVertices, 2] = Double(modelSpaceHPoint.z)
                processedVertices += 1
            }
            
            let meshIndexBuffers = indexBuffers[meshNodeIndex]
            let meshIndexCounts = indexCount[meshNodeIndex]
            
            for i in 0..<meshIndexBuffers.count {
                let indexBuffer = meshIndexBuffers[i]
                let indexCount = meshIndexCounts[i]
                
                assert(indexCount % 3 == 0)
                let triangleCount = indexCount / 3
                faces.resize(triangleCount, 3)
                let indexBufferAddress = indexBuffer.contents().advanced(by: indexBuffersOffsets[meshNodeIndex][i])
                
                // loop over triangles
                for t in 0..<triangleCount {
                    let i: UInt32 = indexBufferAddress.assumingMemoryBound(to: UInt32.self).advanced(by: 3 * t).pointee
                    let j: UInt32 = indexBufferAddress.assumingMemoryBound(to: UInt32.self).advanced(by: 3 * t + 1).pointee
                    let k: UInt32 = indexBufferAddress.assumingMemoryBound(to: UInt32.self).advanced(by: 3 * t + 2).pointee
                    
                    (faces.valuesPtr.pointer + (3 * t)).initialize(to: Int(i))
                    (faces.valuesPtr.pointer + (3 * t + 1)).initialize(to: Int(j))
                    (faces.valuesPtr.pointer + (3 * t + 2)).initialize(to: Int(k))
                    /*
                    (faces.valuesPtr.pointer + (3 * t)).initialize(
                        from: indexBufferAddress.assumingMemoryBound(to: UInt32.self).advanced(by: 3 * t),
                        count: 3)*/
                }
            }
        }
    }
}
