
src = @__DIR__

if !("Bot.txt" in readdir(src))
    touch("$(src)/Bot.txt")
end

f = open("$(src)/Bot.txt")
const CHAT_ID = readline(f)
const BOT_TOKEN = readline(f)
tg = TelegramClient(BOT_TOKEN, chat_id = CHAT_ID)

function sendUpdate(update::String)
    Telegram.sendMessage(tg, text=update)
    nothing
end
