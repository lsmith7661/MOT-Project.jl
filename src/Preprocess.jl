function Preprocess(VideoString::String)
    AV = VideoIO.open(VideoString)
    videoObj = VideoIO.openvideo(AV)
    raw = VideoIO.read(videoObj)
    mytype = typeof(Gray.(raw))
    img = convert(Array{Float64},Gray.(raw)).*255

    ImgSet = [img]
    counter = 0
    while !eof(videoObj)
        read!(videoObj, raw)
        img = convert(Array{Float64},Gray.(raw)).*255
        push!(ImgSet,img)
        counter+=1
    end
    close(videoObj)

    return ImgSet
end

export Preprocess
